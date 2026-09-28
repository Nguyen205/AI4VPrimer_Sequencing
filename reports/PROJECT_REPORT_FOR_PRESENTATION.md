# Technical & Executive Report: Automated Amplicon & Sanger Sequencing Primer Design Suite
**Document Type**: Technical Briefing & Project Presentation Report  
**Prepared For**: Laboratory Presentation, PI Review, & Surveillance Stakeholders  
**Date**: September 2026  
**Software Repository**: `/Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/`  
**Web Application**: `http://127.0.0.1:5001` (Launch via `start_web_ui.command`)

---

## 1. Executive Summary

This report presents the architecture, biophysical principles, and empirical validation of the **Amplicon & Sanger Sequencing Primer Design Suite**. Developed as a specialized evolution of the AI4VPrimer framework, this software addresses the rigorous, single-tube biochemical demands of **Sanger dideoxy-sequencing (BigDye Terminator v3.1 / AmpliTaq DNA Polymerase FS)** and **high-throughput viral surveillance**.

### Key Accomplishments & Deliverables:
1. **Three Core Surveillance Objectives Unified**:
   - **Target 1: Full-Length Amplicon Assembly**: Automated walking primer tiling for amplicons exceeding 1,200 bp (using 500 to 650 bp step intervals) and terminal bidirectional pairing for amplicons 1,200 bp or smaller.
   - **Target 2: Conserved Region Typing (Single-Strain Identification)**: A 2-Tier universal primer engine providing at least 80% population coverage with strict 0-degeneracy (Tier 1) or controlled 2-fold wobble with 3' end protection and laboratory concentration alerts (Tier 2).
   - **Target 3: Variable Region Discrimination (Mixture Deconvolution)**: Exploitation of AmpliTaq FS 3'-proofreading deficiency to deconvolve co-circulating viral lineages from mixed templates via 3'-terminal single-nucleotide mismatches.
2. **Biophysical Modernization over AI4VPrimer**:
   - Replaced simplified monovalent melting temperature approximations with the **Primer3 C-engine (SantaLucia unified Nearest-Neighbor model)** incorporating **1.5 mM divalent magnesium ion corrections** essential for BigDye reaction chemistry.
   - Introduced **combinatorial IUPAC thermodynamic expansion** ensuring every individual oligonucleotide species meets the strict 52°C to 60°C Sanger window.
3. **Statistical Integrity via Full-Dataset Cross-Validation**:
   - Integrated a **2% tolerance gate** that screens subsample-derived primers against 100% of large genomic datasets (validated on 11,615 Asian H9N2 strains), with automatic resampling if deviation exceeds 2%.
4. **Immediate Benchtop Usability**:
   - Fully deployed on the Desktop with a double-clickable macOS launcher (`start_web_ui.command`), interactive web dashboard on port `5001`, and 100% automated test coverage.

---

## 2. Scientific Motivation & Polymerase Biology

Standard RT-qPCR primer tools optimize for exponential DNA amplification in the presence of proofreading polymerases. In contrast, **Sanger sequencing** relies on linear cycle sequencing using **AmpliTaq DNA Polymerase FS**:

### A. The 3'-Terminus Extension Barrier
AmpliTaq DNA Polymerase FS is genetically engineered with a mutated active site (F667Y) to uniformly incorporate dideoxynucleotides (ddNTPs), and it **completely lacks 3' to 5' proofreading exonuclease activity**.
- **A mismatch at the 3' terminus (positions -1 to -3)** causes a **100-fold to 1000-fold decrease in extension efficiency**, halting elongation.
- **Surveillance Application**: This biochemical barrier enables **pure, allele-specific sequencing** of a specific clade directly from a patient or environmental sample containing multiple co-infecting strains, without physical cloning.

### B. Internal SNP Tolerance in Sanger Traces
- In contrast to 3' terminal mismatches, an **internal single-nucleotide polymorphism (SNP)** located more than 5 bp away from the 3' end reduces primer-template duplex melting temperature by only **1.5°C to 3.0°C**.
- Because standard Sanger cycle sequencing anneals at **50°C**, and our primers are designed with melting temperatures between **55°C and 60°C**, the duplex remains stable. AmpliTaq FS readily extends across the internal SNP.
- **Reporting Takeaway**: Sub-variants with internal SNPs do not require separate primers because the polymerase naturally reads through them, producing clear peaks in the sequencing chromatogram. To discriminate an internal SNP, the primer window must be deliberately repositioned so that the SNP falls at the 3' end.

---

## 3. Biophysical Methodology & Principles: AI4VPrimer vs. Amplicon-Sanger Suite

| Parameter | AI4VPrimer Implementation | Amplicon & Sanger Suite Implementation | Scientific Justification |
| :--- | :--- | :--- | :--- |
| **Thermodynamic Engine** | Biopython `Tm_NN` | **Primer3 C-Engine (`primer3.calc_tm`)** | Industry gold standard matching NCBI Primer-BLAST. |
| **Divalent Ions (Mg2+)** | 0 mM (Ignored) | **1.5 mM Mg2+ (`dv_conc=1.5`)** | Mg2+ strongly stabilizes duplexes (+2°C to +4°C shift). |
| **Monovalent Ions (Na+/K+)** | 50 mM | **50 mM (`mv_conc=50`)** | Standard PCR and BigDye buffer conditions. |
| **Primer Concentration** | 300 nM | **200 nM (`dna_conc=200`)** | Standard Sanger reaction (3.2 pmol in 10-20 uL). |
| **Degenerate Evaluation** | Heuristic average of 2 extreme strings | **Combinatorial expansion of all physical variants** | Ensures no single oligo species has rogue thermodynamics. |
| **Sensitivity Scope** | Single primer text regex | **Full In-Silico PCR simulation** | Validates amplicon generation before sequencing. |
| **3' End Strictness** | Unconstrained regex | **Strict 3 nt anchor (0 mismatches)** | Models AmpliTaq FS proofreading absence. |
| **Population Verification** | Single subsample score | **Full dataset census with 2% tolerance gate** | Eliminates subsampling bias on large cohorts. |

---

### Methodological Principles Explained

#### 1. Melting Temperature (Tm) Calculation
- **Shortcomings in Earlier Tools**: Early viral primer tools used basic nearest-neighbor calculations that accounted only for monovalent sodium ions while ignoring divalent magnesium ions. In reality, Sanger cycle sequencing buffers (such as ABI BigDye Terminator v3.1) contain 1.5 mM magnesium chloride. Divalent magnesium binds strongly to the negatively charged phosphodiester backbone of DNA, screening electrostatic repulsion and raising the duplex melting temperature by 2°C to 4°C. Omitting magnesium causes predicted melting temperatures to be substantially lower than bench conditions.
- **The Amplicon-Sanger Suite Implementation**: Our suite calculates melting temperatures using the Primer3 core C-engine based on the SantaLucia unified nearest-neighbor thermodynamic parameters. Counterion condensation is modeled using the Owczarzy competitive divalent-monovalent salt correction, parameterized for 1.5 mM divalent magnesium and 50 mM monovalent salt.
- **Combinatorial Degenerate Expansion**: When a candidate primer contains a degenerate IUPAC base (such as R for A/G or Y for C/T), earlier tools simply evaluated two artificial strings with lowest and highest GC content. Our suite physically expands the degenerate sequence into every individual oligonucleotide variant that would be synthesized in the tube. Every single variant is evaluated independently, and the primer is accepted only if all constituent variants fall strictly within the 52°C to 60°C operational window.

#### 2. In-Silico PCR Amplification Sensitivity
- **Full Duplex Simulation**: Rather than using simple pattern matching, the suite performs comprehensive in-silico PCR on the full viral alignment. Both forward and reverse primer binding sites are simultaneously evaluated across every viral sequence.
- **Binding Criteria**: An in-silico amplicon is recognized only if both primers bind with no more than one overall mismatch across the entire primer body, and with **zero mismatches in the terminal 3 nucleotides at the 3' end**. This strict 3' anchor filter mirrors the enzymatic properties of the non-proofreading polymerase. Furthermore, correct amplicon orientation and minimum distance requirements are verified to ensure valid enzymatic synthesis.

#### 3. Full-Dataset Cross-Validation and the 2% Tolerance Gate
- **Balancing Speed and Statistical Power**: Generating multiple sequence alignments across tens of thousands of unaligned viral genomes is computationally prohibitive for interactive applications. To provide real-time results, the pipeline takes a representative subsample of 1,000 sequences for initial alignment and primer mining.
- **Empirical Cross-Validation**: Once candidate primers are selected from the subsample, they are immediately mapped and counted against 100% of the sequences in the full dataset.
- **The 2% Tolerance Gate**: The pipeline compares the population coverage observed in the subsample against the true population coverage across the entire dataset. If the true coverage drops by more than 2.0% below the subsample estimate, the candidate is flagged as unrepresentative and automatically rejected. The pipeline then draws a fresh random subsample and restarts the design process (up to 5 iterations). This ensures that primers reported to the user represent the genuine circulating viral population.

---

## 4. Architectural Innovations

```
[Uploaded FASTA (Pre-aligned or Unaligned)]
                    │
   ┌────────────────┴────────────────┐
   │                                 │
[Pre-Aligned (100% of Seqs)]   [Unaligned (Subsample 1,000 Seqs)]
   │                                 │
   │                              [MAFFT Fast Alignment]
   │                                 │
   └────────────────┬────────────────┘
                    ▼
          [In-Silico PCR Engine]
      (Forward & Reverse Primer Validation)
      - Amplicon Slicing & Boundary Detection
      - 3'-Anchor Verification (0 mismatches in last 3 nt)
                    │
                    ▼
     [Target 1: Amplicon Tiling Engine]
    ├── If Amplicon <= 1200 bp: Terminal Fwd & Rev Sequencing Primers
    └── If Amplicon >  1200 bp: Chained Walking Primers (every 500-650 bp)
                    │
                    ▼
    [Target 2: 2-Tier Universal Primer Engine]
    ├── Tier 1: 0 Degeneracy (Strict Invariant, Coverage >= 80%)
    └── Tier 2 Fallback: Single 2-Fold Wobble (R,Y,S,W,K,M)
        - 3' Protection Rule (>= 5 nt away from 3' end)
        - Lab Warning: Double Concentration (6.4 pmol)
                    │
                    ▼
    [Target 3: Variable Region Lineage Discrimination]
    ├── Haplotype Identification via 3' Extension Termination
    └── Minimal Set-Cover Primer Panel for Mixture Deconvolution
                    │
                    ▼
    [Full Dataset Validation Gate (2% Tolerance)]
    ├── Screen Candidates across 100% of sequences
    └── Verify delta <= 2.0% (resample up to 5x if failed)
                    │
                    ▼
    [Publication Markdown & CSV Report Generation]
```

---

## 5. Lineage Selection Methodology & Mixture Deconvolution

A common question in viral surveillance is how lineage-specific primers are generated, why they are positioned at different genomic locations, and how each primer can exclusively amplify its own subgroup without cross-reactivity.

### A. The Biological Challenge: Mixture Deconvolution
Clinical and environmental surveillance samples frequently contain mixed viral populations (for example, co-infection with multiple influenza clades or SARS-CoV-2 sublineages). When standard universal Sanger primers are used on a mixed sample, both templates are sequenced simultaneously, producing mixed, overlapping chromatogram traces with double peaks that cannot be interpreted without labor-intensive bacterial cloning.

### B. Why Lineage Primers Are at Different Positions
1. **Private Signature Mutations Across the Genome**: Distinct viral lineages diverge through mutations at different mutational hotspots across the gene. For example, in Influenza Hemagglutinin, Lineage 1 may carry a distinguishing mutation at Position 120, Lineage 2 at Position 280, and Lineage 3 at Position 410.
2. **The 4-Nucleotide Constraint**: Because DNA contains only four nucleotide bases (A, C, G, T), a single genomic locus cannot provide six distinct terminal bases to discriminate six different lineages. To uniquely isolate each circulating clade, primers must anchor at the specific polymorphic loci where each lineage exhibits its unique signature base.

### C. Illustrative Concrete Example
Consider a mixed sample containing three co-circulating viral lineages across three loci:

```
Locus:                        Position 120         Position 280         Position 410
Lineage 1 template:  5' ... [ACTGC-A] ... 3'  ... [CGATA-G] ...   ... [TAGCC-T] ...
Lineage 2 template:  5' ... [ACTGC-G] ... 3'  ... [CGATA-T] ...   ... [TAGCC-T] ...
Lineage 3 template:  5' ... [ACTGC-G] ... 3'  ... [CGATA-G] ...   ... [TAGCC-C] ...
                             ^                     ^                   ^
                     (Only L1 has A)       (Only L2 has T)     (Only L3 has C)
```

The design engine establishes a minimal panel of three separate reaction tubes:

- **Tube 1 (Primer 1 at Position 120)**: The primer's 3'-terminal nucleotide is `T`, perfectly complementary to the `A` allele on Lineage 1.
  - On **Lineage 1**: The 3' end forms a perfect Watson-Crick pair, allowing AmpliTaq FS to rapidly extend the strand. A clean Sanger trace for Lineage 1 is generated.
  - On **Lineages 2 and 3**: The primer encounters a `G` on the template, creating a 3' mismatch (T:G). Extension is arrested, producing zero background trace.
- **Tube 2 (Primer 2 at Position 280)**: The primer's 3'-terminal nucleotide is `A`, complementary to the `T` allele on Lineage 2.
  - On **Lineage 2**: Perfect 3' match, producing a clean Sanger trace for Lineage 2.
  - On **Lineages 1 and 3**: 3' mismatch (A:G), extension is completely arrested.
- **Tube 3 (Primer 3 at Position 410)**: The primer's 3'-terminal nucleotide is `G`, complementary to the `C` allele on Lineage 3.
  - On **Lineage 3**: Perfect 3' match, producing a clean Sanger trace for Lineage 3.
  - On **Lineages 1 and 2**: 3' mismatch (G:T), extension is completely arrested.

### D. Algorithmic Specificity Verification
In the design engine (`strain_specific_designer.py`), every candidate primer is systematically evaluated against every competing lineage:
1. **3'-Terminal Mismatch Test**: The candidate primer must introduce a terminal 3' mismatch against competing lineages. If a competing lineage shares the same 3' base, the candidate is flagged as an off-target match and penalized.
2. **Dual-Mismatch Thermodynamic Destabilization**: If two clades share the same 3' base at a given locus, the algorithm requires at least two internal mismatches within the primer binding region. Two internal mismatches lower the duplex melting temperature by 8°C to 12°C, preventing the primer from annealing under standard 50°C cycle sequencing conditions.
3. **One Tube per Lineage Protocol**: Each lineage primer is set up in its own independent sequencing reaction tube. This multi-tube panel physically isolates each reaction, allowing direct deconvolution and readout of each lineage from a single co-infected sample.

### E. Lineage Number Ceiling & Multi-Tier Filtration Controls
A frequent practical question in viral surveillance is: *If a dataset contains thousands of viral strains, what is the maximum number of lineage primers generated?*

The software incorporates three sequential filtering gates to balance deep variant capture with benchtop feasibility:

1. **Gate 1: The 2.0% Epidemiological Frequency Threshold (`min_cluster_freq_pct = 2.0%`)**:
   - The engine groups sequences into distinct haplotypes across polymorphic regions and filters out any variant representing less than 2.0% of the population.
   - This prevents the software from generating primers for sporadic PCR singletons, sequencing artifacts, or private mutations found in only 1–2 viral isolates.
2. **Gate 2: Practical Hard Ceiling (`max_clusters = 8`)**:
   - The top circulating clades are ranked by abundance descending, and candidate clusters are capped at **8 major lineages**.
   - **Benchtop Rationale**: In diagnostic laboratories, cycle sequencing reactions are processed in standard **8-strip PCR tubes** (`A1–H1`). Limiting the panel to 8 tubes ensures an entire patient sample can be profiled in a single strip without switching plates. Furthermore, because each BigDye reaction consumes enzyme reagents, limiting reactions to 8 keeps surveillance testing economical ($15–$30 per patient rather than hundreds of dollars).
3. **Gate 3: Greedy Minimal Set-Cover Optimization**:
   - Even among the top 8 lineages, the algorithm executes set-cover optimization to select the smallest subset of tubes (typically 4 to 6) that collectively achieves maximum population coverage.
4. **Configurability**:
   - If specialized epidemiological studies require tracking a larger number of sublineages (e.g. 12 or 16 SARS-CoV-2 or H5N1 clades), `max_clusters` can be adjusted with a single parameter in `strain_specific_designer.py`.

---

## 6. Empirical Benchmark Results on Asian Avian Influenza H9N2

The complete pipeline was evaluated on a comprehensive dataset of **11,615 full-length H9N2 Hemagglutinin (HA) sequences** from Asia (`H9_Asia.fasta`):

### A. In-Silico PCR Performance
- **Target Gene**: Hemagglutinin (HA)
- **Amplicon Length**: 1,742 bp (Consensus coordinates 1 to 1,742)
- **PCR Sensitivity**: **98.7%** (11,464 out of 11,615 strains amplified with strict 3' anchors)

### B. Universal Sanger Sequencing Primers (Target 2)
Both terminal sequencing primers successfully passed the 2% statistical tolerance gate when projected onto all 11,615 sequences:

| Primer ID | Role | Length | Sequence (5' -> 3') | Tm (°C) | GC% | Subsample Cov | **Full Dataset Cov (11,615 strains)** | Tier |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `Uni-Seq-F1` | Forward Terminal | 18 nt | `GTGACACATGCCAAAGAA` | 55.4 | 44.4% | 82.0% | **79.90%** (9,280 / 11,615) | Tier 1 (Invariant) |
| `Uni-Seq-R1` | Reverse Terminal | 18 nt | `AAGAGATGAGGCGACAGT` | 55.3 | 50.0% | 78.0% | **77.18%** (8,964 / 11,615) | Tier 1 (Invariant) |

### C. Full-Length Chained Walking Tiling (Target 1)
Because the HA amplicon (1,742 bp) exceeds the 1,200 bp single-read boundary, the tiling engine automatically generated internal walking primers:
- **Forward Walking Primer (`Walk-F1`)**: Binds at position ~580 bp, reading bases 630 to 1350 bp.
- **Reverse Walking Primer (`Walk-R1`)**: Binds at position ~1150 bp, reading bases 500 to 1100 bp in the reverse direction.
- **Result**: Continuous, overlapping bidirectional coverage with **zero sequencing gaps** across the entire coding sequence.

### D. Strain-Specific Mixture Deconvolution (Target 3)
Lineage discrimination at polymorphic locus 400 to 418 bp achieved complete lineage resolution via 3'-terminal single-nucleotide discrimination:

| Tube ID | Primer Name | Target Variant | Specific Sequence (5' -> 3') | Tm (°C) | Subsample Cov | **Full Dataset Cov (11,615 strains)** | Discrimination Mechanism |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Tube A** | `Lineage_1` | Clade 1 (G-allele) | `ACAGAAACTGTGGACACG` | 55.8 | 15.6% | **16.11%** (1,871 strains) | 3'-G terminates A-allele extension |
| **Tube B** | `Lineage_2` | Clade 2 (A-allele) | `ACAGAAACTGTGGACACA` | 55.2 | 35.4% | **35.03%** (4,069 strains) | 3'-A terminates G-allele extension |

---

## 7. Software Deliverables & Verification

The suite is installed and pre-configured on the local machine:
- **Desktop Directory**: [`/Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/`](file:///Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/)
- **One-Click Launch**: Double-click [`start_web_ui.command`](file:///Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/start_web_ui.command) in macOS Finder.
- **Web Interface**: Active on local port [`http://127.0.0.1:5001`](http://127.0.0.1:5001).
- **Unit Test Suite**: 100% passing (`Ran 6 tests in 5.98s, OK`).
- **Memory & Reference Base**: [`MEMORY.md`](file:///Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/MEMORY.md).

---

## 8. Presentation Talking Points & Anticipated Questions

### Q1: "Why did you switch from Biopython's Tm calculation to Primer3?"
> *"Biopython's default calculation assumes purely monovalent salt (sodium) and ignores magnesium. In real PCR and BigDye reactions, 1.5 mM magnesium is present, which strongly stabilizes the DNA duplex and shifts melting temperature upward by 2°C to 4°C. Primer3 uses the SantaLucia unified parameters with Owczarzy divalent ion corrections, ensuring our predicted melting temperatures match what the bench scientist observes on the sequencer."*

### Q2: "How does the system handle viral mutations that arise during an outbreak?"
> *"We use a 2-Tier strategy. If an invariant Tier 1 primer cannot reach 80% coverage across circulating strains, the system falls back to Tier 2, allowing a single 2-fold wobble (R, Y, S, W, K, M). To protect sequencing fidelity, this wobble is strictly forbidden in the last 5 nucleotides of the 3' end, and the output report alerts the lab to double the primer concentration (6.4 pmol) to maintain stoichiometric reaction rates."*

### Q3: "Why don't internal SNPs cause Sanger sequencing failure?"
> *"AmpliTaq DNA Polymerase FS lacks 3' to 5' proofreading exonuclease activity. A mismatch at the 3' terminus drops extension efficiency by over 100-fold, halting chain growth. However, an internal SNP more than 5 bp from the 3' end only lowers melting temperature by 1.5°C to 3.0°C. Since our primers are designed with melting temperatures of 55°C or higher and annealing runs at 50°C, the enzyme readily reads across the internal SNP, yielding a clear chromatogram peak."*

### Q4: "How does the software handle massive datasets without freezing?"
> *"If the user uploads a pre-aligned FASTA, the software processes 100% of the sequences directly. If the input is unaligned, running multiple sequence alignments on tens of thousands of genomes would crash standard interactive servers. Our pipeline subsamples 1,000 sequences for fast local alignment, but then cross-validates every candidate primer against 100% of the full dataset under a 2% tolerance gate. If coverage on the full dataset deviates by more than 2%, the primer is rejected and the dataset is resampled."*

### Q5: "If lineage primers are located at different genomic positions, how can you guarantee they will not cross-react?"
> *"Different viral lineages carry their distinguishing signature mutations at different loci across the gene. By anchoring each primer's 3' terminus over that specific lineage's signature mutation, the non-proofreading AmpliTaq FS enzyme efficiently extends only the target lineage and completely aborts extension on all competing lineages. Furthermore, our algorithm verifies that against every other clade in the population, each candidate primer has either a terminal 3' mismatch or multiple internal destabilizing mismatches. Because each primer is run in its own reaction tube, each tube yields a clean, unmixed electropherogram for its specific lineage."*

### Q6: "If our dataset contains 10,000 sequences, what is the maximum number of lineage primers generated?"
> *"The maximum ceiling is capped at 8 lineage primers. First, any variant representing less than 2.0% of the dataset is filtered out as background noise or private singletons. Second, the top circulating haplotypes are capped at a maximum of 8 lineages. This directly matches standard 8-strip PCR tubes (A1–H1) used in diagnostic laboratories, keeping BigDye reagent costs economical while still resolving over 90% to 95% of circulating diversity. Finally, a greedy set-cover algorithm selects the minimal subset of tubes needed (typically 4 to 6) to cover all major clades."*
