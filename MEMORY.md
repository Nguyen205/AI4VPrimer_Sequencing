# Amplicon & Sanger Sequencing Primer Design Suite: Project Memory & Knowledge Base
**Last Updated**: September 2026  
**Project Path**: `/Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/`  
**Dev Repository**: `/Users/aijiazhou/AmpliconSangerPrimer/`  
**Web UI**: `http://127.0.0.1:5001` (Launch via `start_web_ui.command` or `python3 app.py`)

---

## 1. Project Background & Motivation
This project develops an automated, enterprise-grade **Amplicon & Sanger Sequencing Primer Design Suite** for viral surveillance, inspired by and expanding upon [Nguyen205/AI4VPrimer_Influenza](https://github.com/Nguyen205/AI4VPrimer_Influenza).

While AI4VPrimer focused primarily on standard RT-PCR amplification, Sanger sequencing places much stricter biophysical and biochemical constraints on primer design:
- Requires single-tube dideoxy-terminator chemistry (BigDye Terminator v3.1 / AmpliTaq DNA Polymerase FS).
- High sensitivity to primer-dimers, secondary hairpins, and template secondary structures.
- Requires high sequence conservation across viral clades while enabling targeted discrimination of lineages and mixtures.
- Demands tailored walking strategies depending on amplicon length ($\le 1200\text{ bp}$ vs. $> 1200\text{ bp}$).

---

## 2. Chronological Conversation & Milestone Log

### Milestone 1: Core Architecture & 3 Target Objectives
- **Target 1: Full-Length Amplicon Assembly**: Tiling and walking primers to sequence complete viral genes (e.g. Influenza HA, NA) without gaps.
- **Target 2: Conserved Region Typing (Single-Strain Identification)**: High-conservation universal primers capturing maximal strains across historical and emerging clades.
- **Target 3: Variable Region Discrimination (Mixture Deconvolution)**: Strain- and lineage-specific primers exploiting terminal single-nucleotide differences to isolate specific variants within co-infected mixtures.

### Milestone 2: Biophysical Standards for Sanger Primers
- **Primer Length**: $18 - 24\text{ nt}$ (target $20\text{ nt}$).
- **Melting Temperature ($T_m$)**: Strict $52^\circ\text{C} - 60^\circ\text{C}$ (optimum $55^\circ\text{C} - 58^\circ\text{C}$). Paired walking primers must have $\Delta T_m \le 3^\circ\text{C}$.
- **GC Content**: $40\% - 60\%$.
- **Homopolymer Penalty**: Reject any candidate with $\ge 4$ consecutive identical bases (`AAAA`, `GGGG`, etc.) to prevent polymerase slippage.
- **Secondary Structures**: Hairpin $T_m < 45^\circ\text{C}$, homodimer $T_m < 40^\circ\text{C}$, heterodimer $T_m < 40^\circ\text{C}$ (evaluated via `primer3-py`).

### Milestone 3: 2-Tier Universal Primer Strategy
1. **Tier 1 (Strict Invariant)**: 0 degenerate bases. Requires $100\%$ exact match across target sequences. Preferred for clean Sanger traces.
2. **Tier 2 (Controlled Wobble Fallback)**:
   - Activated automatically only when Tier 1 coverage drops below $80\%$.
   - Allows **at most one 2-fold degenerate IUPAC base** (`R, Y, S, W, K, M`).
   - **3'-Protection Rule**: Degenerate base must be $\ge 5\text{ bp}$ away from the 3' terminus to safeguard primer extension efficiency.
   - **Laboratory Alert**: Automatically flags the primer report:
     > `DOUBLE PRIMER CONCENTRATION (6.4 pmol) IN SANGER MIX`: When using a 2-fold degenerate primer in BigDye chemistry, standard primer amount (3.2 pmol) must be doubled to 6.4 pmol so that each individual oligonucleotide variant reaches full working concentration.
   - **`primer3-py` IUPAC Handling**: Since `primer3-py` rejects IUPAC codes, `expand_iupac_variants()` expands the degenerate primer into all constituent oligonucleotides, computing thermodynamics across each variant and selecting the worst-case thermodynamic penalty for safety.

### Milestone 4: Amplicon Tiling & Walking Logic
- **Amplicons $\le 1200\text{ bp}$**:
  - Sanger read length is typically $700 - 900\text{ bp}$.
  - Forward and Reverse terminal sequencing primers (nested $\sim 50\text{ bp}$ inside PCR primers) provide overlapping bidirectional coverage of the entire amplicon without requiring internal walking primers.
- **Amplicons $> 1200\text{ bp}$**:
  - Automatically tiles internal walking primers spaced every **$500 - 650\text{ bp}$**.
  - Includes universal M13 sequencing tails (Forward: `TGTAAAACGACGGCCAGT`, Reverse: `CAGGAAACAGCTATGACC`) option for standardized high-throughput core facility sequencing.

### Milestone 5: Pre-Aligned vs. Unaligned FASTA Handling
- **Pre-Aligned FASTA (Recommended)**:
  - If user uploads a pre-aligned alignment, the pipeline processes **100% of all sequences** directly. No subsampling is performed.
  - The Web UI displays a prominent banner recommending pre-aligned FASTA uploads for maximum statistical power.
- **Unaligned FASTA Subsampling**:
  - Running MAFFT on $>10,000$ unaligned viral genomes is computationally prohibitive for interactive web requests.
  - **Subsample Size**: Fixed to **1,000 sequences** (benchmarking confirmed the runtime difference between 500 and 1,000 sequences is only a few seconds, while 1,000 significantly improves representation of low-frequency clades).
  - **Full Dataset Validation & 2% Tolerance Gate**: Primers designed on the 1,000-sequence subsample are immediately screened against 100% of the sequences in the full dataset (e.g. all 11,615 strains in `H9_Asia.fasta`). If coverage on the full dataset deviates by $> 2\%$ below the subsample estimate, the pipeline rejects the candidate and resamples (up to 5 attempts).

### Milestone 6: Mixture Deconvolution & Sanger Polymerase Biology
The user raised critical biological questions regarding SNP tolerance and lineage classification:
1. **Does Sanger sequencing require a strict match? How do SNPs affect it?**
   - **AmpliTaq DNA Polymerase FS** (the proprietary Taq mutant in ABI BigDye Terminator kits) lacks $3' \rightarrow 5'$ proofreading exonuclease activity.
   - **3'-Terminus Strictness**: A single mismatch at the 3' terminus (the final $1 - 3\text{ nt}$) causes a **$100\times - 1000\times$ decrease in extension efficiency**, effectively arresting the reaction. This biochemical property is the foundation for our allele-specific lineage primers.
   - **Internal SNP Tolerance**: An internal mismatch ($> 5\text{ bp}$ away from the 3' end) only lowers melting temperature by $\sim 1.5^\circ\text{C} - 3^\circ\text{C}$. Because standard Sanger annealing cycles run at $50^\circ\text{C}$ and our primers have $T_m \ge 55^\circ\text{C}$, the enzyme easily binds and extends across internal SNPs. Thus, internal SNPs are read as clear secondary peaks (or single peaks if homozygous) in the chromatogram.
2. **Why are sub-variants with 1 internal SNP not assigned a new lineage primer?**
   - If an internal SNP does not prevent extension, designing a separate primer for it would result in competitive binding and cross-priming, ruining mixture deconvolution.
   - To discriminate an internal SNP, the primer design window must be shifted so that the SNP is positioned directly at the **3' terminus**.

### Milestone 7: Lineage Classification Census Fix
- **Issue Identified**: In initial testing on `H9_Asia.fasta`, the lineage report only listed 12 and 6 sequences for distinct lineages despite 500 sequences being sampled.
- **Root Cause**: The clustering algorithm had been clustering sequences over a broad $400\text{ bp}$ amplicon window. Natural divergence across those 400 bases fragmented the 500 sequences into dozens of micro-clusters.
- **Resolution**: Updated `pipeline.py` to calculate the **true binding counts** for each specific designed lineage primer at its exact $18 - 24\text{ nt}$ locus across both the subsample and the full dataset.
- **Verified Benchmark (`H9_Asia.fasta`, 11,615 total sequences)**:
  - `Uni-Seq-F1` (`GTGACACATGCCAAAGAA`): 410 / 500 ($82.0\%$) subsample $\rightarrow$ **9,280 / 11,615 ($79.90\%$) full dataset**.
  - `Uni-Seq-R1` (`AAGAGATGAGGCGACAGT`): 390 / 500 ($78.0\%$) subsample $\rightarrow$ **8,964 / 11,615 ($77.18\%$) full dataset**.
  - `Lineage_1` (`ACAGAAACTGTGGACACG`): 78 / 500 ($15.6\%$) subsample $\rightarrow$ **1,871 / 11,615 ($16.11\%$) full dataset**.
  - `Lineage_2` (`ACAGAAACTGTGGACACA`): 177 / 500 ($35.4\%$) subsample $\rightarrow$ **4,069 / 11,615 ($35.03\%$) full dataset**.

### Milestone 8: Lineage Selection Logic Across Different Genomic Loci & Mixture Deconvolution
- **Why Lineage Primers Are at Different Positions**:
  - Different viral lineages mutate and diverge at different mutational hotspots across the gene (private signature mutations).
  - DNA has only 4 bases ($A, C, G, T$). A single locus cannot provide 6 distinct bases to separate 6 lineages. Therefore, primers must anchor at different polymorphic loci across the amplicon where each lineage has its distinct signature base.
- **Biochemical Specificity & Proof of Non-Crossing**:
  - Cycle sequencing utilizes **AmpliTaq DNA Polymerase FS**, which lacks $3' \rightarrow 5'$ proofreading exonuclease activity.
  - A mismatch at the 3' terminus drops extension efficiency by $100\times - 1000\times$, completely arresting extension on non-target lineages.
  - Candidate primers are cross-validated against every other cluster in the dataset. If any competing lineage shares the same 3' base without $\ge 2$ internal destabilizing mismatches ($\Delta T_m \approx 8-12^\circ\text{C}$ penalty), the candidate is rejected.
- **Multi-Tube Workflow**: Each lineage primer is set up in a separate reaction tube (`Tube 1`, `Tube 2`, ..., `Tube 6`). Each tube selectively sequences only its target lineage from a mixed sample, eliminating double chromatogram peaks.
- **Lineage Number Ceiling & Multi-Tier Filtration Logic**:
  - **Gate 1 (2% Frequency Filter)**: Excludes singletons, private mutations, and PCR artifacts (< 2% prevalence).
  - **Gate 2 (Default Ceiling: 8 Lineages)**: Capped at `max_clusters = 8` to perfectly match standard 8-strip PCR tubes (`A1–H1`) and prevent excessive BigDye reagent consumption.
  - **Gate 3 (Greedy Minimal Set Cover)**: Mathematically chooses the minimum tubes (typically 4–6) needed to cover >95% of circulating clades.
  - **Configurability**: `max_clusters` and `min_cluster_freq_pct` can be customized in `strain_specific_designer.py`.
- **Reporting & Presentation Update**:
  - Removed all raw LaTeX math equations from the executive report, replacing them with clear narrative paragraphs explaining nearest-neighbor thermodynamics, $1.5\text{ mM } Mg^{2+}$ salt corrections, in-silico PCR sensitivity, and the 2% tolerance gate.
  - Consolidated redundant position columns into a unified `Position (bp)` column (`Cols X – Y`).
  - Added dedicated Section 5.E and talking point Q6 on lineage number control to the presentation report.
  - Formatted and generated both Markdown (`.md`), Microsoft Word (`.docx`), and Google Docs-ready (`.html`) report formats.

---

## 3. Directory Layout & Module Responsibilities

```
/Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/
├── app.py                          # Flask Web Application (GUI on http://127.0.0.1:5001)
├── start_web_ui.command             # Double-clickable macOS Finder script to launch web UI
├── README.md                       # Comprehensive User Guide & Methodology
├── MEMORY.md                       # This Persistent System Knowledge Base
├── requirements.txt                # Biopython, primer3-py, Flask
├── templates/
│   └── index.html                  # Bootstrap 5 GUI with Pre-Aligned recommendation banner
├── reports/
│   └── sanger_primer_H9_Asia.md    # Publication-ready validation report on 11,615 H9 strains
├── tests/
│   └── test_pipeline.py            # Automated test suite (6/6 unit tests passing)
└── sanger_designer/                # Core Algorithmic Engine
    ├── __init__.py                 # Package exports
    ├── msa_utils.py                # MSA loading, IUPAC processing, 1000-subsample & MAFFT fallback
    ├── insilico_pcr.py             # In-silico PCR amplification, mismatch tolerance, amplicon slicing
    ├── universal_designer.py       # 2-Tier universal primer mining, IUPAC thermodynamic expansion
    ├── strain_specific_designer.py # Haplotype clustering & 3'-terminal extension termination
    ├── tiling.py                   # <=1200 bp terminal vs. >1200 bp walking tiling logic
    └── pipeline.py                 # Master orchestrator, 2% tolerance gate, true census reporting
```

---

## 4. Key Algorithmic Implementations

### A. 2-Tier Universal Primer Discovery (`universal_designer.py`)
- Slides window of length $k \in [18, 24]$ across the consensus sequence.
- Evaluates candidate at position $i$:
  - Calculates exact conservation ratio $C = \frac{\text{Count}(\text{dominant } k\text{-mer})}{\text{Total valid sequences}}$.
  - If $C \ge 0.80$, tests Primer3 biophysical filters ($T_m \in [52, 60]$, GC $\in [40, 60]$, homopolymer $\le 3$, hairpin $T_m < 45$, dimer $T_m < 40$). If passing, marks as **Tier 1 (Strict Invariant)**.
  - If $C < 0.80$, scans for positions with a single binary SNP. If inserting a 2-fold degenerate IUPAC base at position $p$ with $(k - 1 - p) \ge 5$ brings coverage $\ge 0.80$, performs `expand_iupac_variants()`. If all variants pass biophysics, marks as **Tier 2 (Controlled Wobble Fallback)** and attaches the laboratory double-concentration alert.

### B. Mixture Deconvolution via 3' Mismatch (`strain_specific_designer.py`)
- To distinguish co-circulating variants in an amplicon, identifies highly polymorphic loci within the target region.
- Groups sequences by the exact nucleotide at the 3' terminus ($pos = k - 1$).
- Generates allele-specific primers where the $3'$ base is perfectly complementary to Variant A but mismatched against Variant B.
- Because AmpliTaq FS lacks 3' proofreading, Variant B cannot be extended, yielding a pure electropherogram for Variant A even in an unpurified viral mixture.

### C. Full-Dataset Validation & 2% Tolerance Gate (`pipeline.py`)
```python
# Validation on full alignment
subsample_cov = primer_result["subsample_coverage"]
full_cov = calculate_exact_binding_count(primer_seq, full_alignment) / len(full_alignment)

if (subsample_cov - full_cov) > 0.02:
    # Coverage dropped by more than 2% in the full population
    resample_subsample_and_redesign(attempt_max=5)
```

---

## 5. Verification & Testing Checklist
- [x] Pre-aligned FASTA runs on 100% of sequences without subsampling.
- [x] Unaligned FASTA subsamples 1,000 sequences and cross-validates against 100% of dataset.
- [x] 2% tolerance gate triggers automatic resampling if subsample is non-representative.
- [x] $T_m$ strictly confined to $52^\circ\text{C} - 60^\circ\text{C}$ with max homopolymer $\le 3$.
- [x] Tier 2 degenerate primers alert laboratory to double primer concentration (6.4 pmol).
- [x] Primer3 IUPAC expansion correctly validates thermodynamics for degenerate mixtures.
- [x] True sequence census reported for all lineages across subsample and full dataset.
- [x] Double-clickable `start_web_ui.command` launches Flask app on port 5001.
- [x] Unit test suite passes 100% (`pytest tests/test_pipeline.py`).
- [x] Alignment columns and amplicon positions explicitly displayed across all report tables (Tiling, Universal Primers, and Lineage Specific Panels).
- [x] Dynamic report tables: Pre-aligned FASTA (100% evaluated) automatically omits duplicate 'Subsample' headers and displays clean 'Population Share' and 'Sequence Count'.
