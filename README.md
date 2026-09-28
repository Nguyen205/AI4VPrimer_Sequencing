# AI4VPrimer Amplicon-Sanger Suite

Automated Sanger sequencing primer design for viral surveillance amplicons, strain identification, and mixture deconvolution.

---

## Key Features

1. **Two-Tier Universal Primer Engine (Single-Strain Testing)**
   - **Tier 1 (Strict Invariant, 0 Degeneracy)**: Identifies 100% conserved sites across input strains with optimal Tm ($55^\circ\text{C} - 60^\circ\text{C}$) and absence of homopolymer repeats ($\ge 4$ identical nt).
   - **Tier 2 (Fallback: 2-Fold Degeneracy)**: If strain diversity limits Tier 1 coverage below $80\%$, automatically enables a single 2-fold IUPAC wobble (`R, Y, S, W, K, M`).
   - **3'-End Protection Rule**: Zero wobble allowed within the last 5 nucleotides of the 3' end to prevent polymerase stalling.
   - **Concentration Warning**: Automatically flags Tier 2 primers with instructions to **double primer concentration (6.4 pmol)** to compensate for 50% signal dilution.

2. **Strain-Specific Discriminatory Panel (Mixture Deconvolution)**
   - When samples are suspected mixtures of strains (co-infections or quasispecies), universal primers cause overlapping chromatogram peaks.
   - Mines allele-specific primers where the **3'-terminal base mismatches off-target clades**, terminating unwanted extension.
   - **Lineage Number & Multi-Tier Filtration Controls**:
     - **Gate 1: 2.0% Epidemiological Frequency Filter**: Excludes rare private mutations, sequencing artifacts, and singleton strains (< 2% prevalence).
     - **Gate 2: Practical Ceiling (Default: 8 Lineages)**: Automatically caps clusters at 8 major lineages, perfectly matching standard laboratory **8-strip PCR tubes** (`A1–H1`) and preventing expensive reagent waste.
     - **Gate 3: Greedy Minimal Set-Cover**: Mathematically selects the minimal number of reaction tubes (typically 4 to 6) needed to achieve maximal population coverage.
     - **Configurable**: Both `max_clusters` (default: 8) and `min_cluster_freq_pct` (default: 2.0%) can be adjusted in `strain_specific_designer.py`.

3. **Amplicon Sizing & Tiling Control**
   - **$\le 1200\text{ bp}$**: Single-step terminal reactions suffice (opposing Forward and Reverse reads overlap seamlessly).
   - **$> 1200\text{ bp}$**: Automatically activates internal primer walking ($\sim 500-650\text{ bp}$ stepping) to generate chained sequencing tiles across long amplicons.

4. **Integrated In-Silico PCR**
   - Evaluates user-provided external PCR primers or auto-slices aligned FASTA amplicons.
   - Calculates true PCR sensitivity across all target strains.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/Pluto-Aijia/AI4VPrimer_Sequencing.git
cd AI4VPrimer_Sequencing

# Install dependencies
pip install -r requirements.txt
```

### Dependencies
- Python 3.8+
- `biopython >= 1.79`
- `numpy >= 1.21`
- `primer3-py >= 0.6.1`
- `flask >= 3.0`

---

## Running the Web Interface

```bash
python app.py
```
Open **`http://127.0.0.1:5001`** in your browser. Configure your aligned FASTA path, set PCR primers, adjust thresholds, and click **Run Sanger Design Engine**.

---

## Python API Usage

```python
from sanger_designer.pipeline import SangerAmpliconPipeline

pipeline = SangerAmpliconPipeline(
    fasta_path="/path/to/aligned_influenza.fasta",
    fwd_pcr_primer="ATGCGATCGATCGATCGATCGA",
    rev_pcr_primer="CGTAACGTAGCTAGCTAGCTA",
    min_coverage_pct=80.0,
    min_tm=55.0,
    max_tm=60.0,
    output_report_path="sanger_primer_report.md"
)

results = pipeline.run()
print(results["report_content"])
```

---

## Directory Structure

```
Amplicon_Sanger_Primer_Suite/
├── app.py                      # Flask web server & API
├── start_web_ui.command        # One-click desktop launcher (macOS)
├── requirements.txt            # Python dependencies
├── README.md                   # Documentation
├── templates/
│   └── index.html              # Interactive browser interface
├── sanger_designer/
│   ├── __init__.py             # Package exports
│   ├── msa_utils.py            # Alignment parsing, coordinate mapping & IUPAC logic
│   ├── insilico_pcr.py         # In-silico PCR amplification & sensitivity calculation
│   ├── universal_designer.py   # Tier 1 & Tier 2 universal Sanger primer design
│   ├── strain_specific_designer.py # Allele-specific 3' mismatch mixture panel
│   ├── tiling.py               # <= 1200 bp single vs > 1200 bp walking logic
│   └── pipeline.py             # Orchestrator & Markdown report generator
└── tests/
    ├── data/                   # Test datasets and synthetic alignments
    └── test_pipeline.py        # Comprehensive test suite
```
