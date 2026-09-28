# Amplicon & Sanger Sequencing Primer Design Report

Generated automatically by **AI4VPrimer Amplicon-Sanger Engine**.

---

## 1. PCR Amplicon Summary
- **Mode:** Whole-Gene Alignment / Auto Mode (Path B)
- **Template FASTA:** `/Users/aijiazhou/Desktop/Primer_Vietnam/H9_Asia.fasta` (500 sequences)
- **Forward PCR Primer:** `Auto / Full-Length`
- **Reverse PCR Primer:** `Auto / Full-Length`
- **In-Silico PCR Sensitivity:** **100.0%** (500/500 sequences amplified)
- **Amplicon Span:** MSA columns `0` to `1759`
- **Amplicon Length:** **1759 bp**

---

## 2. Amplicon Sizing & Tiling Strategy
- **Amplicon Length:** 1759 bp
- **Tiling Required:** `YES (> 1200 bp)`
- **Strategy:** Multi-Primer Walking (4 primers spaced every ~550 bp)
- **Details:** Amplicon length (1759 bp) exceeds the 1200 bp threshold. 4 chained sequencing primers designed to achieve full-length coverage.

### Tiled Primer Walking Schedule
| Step | Direction | Start (bp) | End (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Subsample Cov | Full Dataset Cov | Tier | Covered Window |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | FORWARD | 115 | 135 | 20 nt | `GGCTAYCAATCAACAAACTC` | 55.34 | 40.0% | 80.8 % | **78.61 %** | Tier 2 (2-fold) | 115 - 815 bp |
| 3 | FORWARD | 452 | 471 | 19 nt | `CAGACACRATCTGGAATGT` | 55.61 | 42.11% | 77.0 % | **73.13 %** | Tier 2 (2-fold) | 452 - 1152 bp |
| 4 | FORWARD | 972 | 991 | 19 nt | `ATGCATTTGGAAACTGCTC` | 56.27 | 42.11% | 72.8 % | **69.36 %** | Tier 1 (0 degen) | 972 - 1672 bp |
| 4 | REVERSE | 1634 | 1652 | 18 nt | `AAGAGATGAGGCGACAGT` | 57.35 | 50.0% | 78.0 % | **77.18 %** | Tier 1 (0 degen) | 952 - 1652 bp (Reverse) |

---

## 3. Mode 1: Universal Sanger Sequencing Primers (Single-Strain Testing)
Use these primers first. If your sample contains a single pure strain, these will yield clean, sharp chromatograms.

### Forward Sequencing Primer (`Uni-Seq-F1`)
- **Tier:** 🟢 Tier 1 (Strict Invariant, 0 Degeneracy)
- **Sequence (5' → 3'):** `GTGACACATGCCAAAGAA`
- **Length:** 18 nt | **Tm:** 55.26 °C | **GC:** 44.44%
- **Discovery Subsample Coverage:** 82.0%
- **FULL DATASET GROUND TRUTH COVERAGE:** **79.9%** (9280 / 11615 strains in `H9_Asia.fasta`)
- **Secondary Structures:** Hairpin $\Delta G$: 0.0 kcal/mol | Dimer $\Delta G$: -4.11 kcal/mol
- **Concentration Protocol:** Standard 1x (e.g. 3.2 pmol per 10 uL reaction)

### Reverse Sequencing Primer (`Uni-Seq-R1`)
- **Tier:** 🟢 Tier 1 (Strict Invariant, 0 Degeneracy)
- **Sequence (5' → 3'):** `AAGAGATGAGGCGACAGT`
- **Length:** 18 nt | **Tm:** 57.35 °C | **GC:** 50.0%
- **Discovery Subsample Coverage:** 78.0%
- **FULL DATASET GROUND TRUTH COVERAGE:** **77.18%** (8964 / 11615 strains in `H9_Asia.fasta`)
- **Secondary Structures:** Hairpin $\Delta G$: 0.0 kcal/mol | Dimer $\Delta G$: -3.9 kcal/mol
- **Concentration Protocol:** Standard 1x (e.g. 3.2 pmol per 10 uL reaction)

---

## 4. Mode 2: Strain-Specific Panel (Mixture Deconvolution)
If your clinical/field sample is a **mixture of strains**, universal primers will produce double peaks. 
Use the minimal panel below. Each primer selectively amplifies its specific clade via **3'-terminal extension termination**.

### How Lineages are Determined
1. **Target Region Partitioning:** The alignment is clustered by shared nucleotide motifs across the variable target domain, grouping identical alleles into biological lineages.
2. **3'-Terminal Selectivity:** Primers are placed so the very last 3' nucleotide matches only the intended lineage with 0 mismatches, but mismatches all other lineages, arresting non-specific Taq polymerase extension.

- **Clusters/Variants Discovered:** 2
- **Panel Tubes Required:** 2
- **Total Cumulative Coverage:** **51.0%** (Subsample) | **51.14%** (Full Dataset)

### Lineage Breakdown in Full Dataset
| Lineage ID | Subsample Count | Subsample Share | Full Dataset Count | Full Dataset Share | Assigned Test Tube |
|---|---|---|---|---|---|
| **Lineage_1** | 78 / 500 seqs | 15.6% | 1871 / 11615 seqs | **16.11%** | **Tube_1** (`Spec_F_Lineage_1`) |
| **Lineage_2** | 177 / 500 seqs | 35.4% | 4069 / 11615 seqs | **35.03%** | **Tube_2** (`Spec_F_Lineage_2`) |

### Minimal Panel Tubes
| Tube ID | Primer ID | Target Variant | Sequence (5' → 3') | Length | Tm (°C) | Subsample Share | Full Dataset Share | Reaction Protocol |
|---|---|---|---|---|---|---|---|---|
| **Tube_1** | `Spec_F_Lineage_1` | Lineage_1 | `ACAGAAACTGTGGACACG` | 18 nt | 56.95 | 15.6% | **16.11%** | Test in Tube 1. Selective 3' extension ensures clean trace for Lineage_1 even in mixtures. |
| **Tube_2** | `Spec_F_Lineage_2` | Lineage_2 | `ACAGAAACTGTGGACACA` | 18 nt | 55.64 | 35.4% | **35.03%** | Test in Tube 2. Selective 3' extension ensures clean trace for Lineage_2 even in mixtures. |

---

## 5. Oligo Synthesis Order Sheet
Copy and paste these sequences directly into your synthesis provider (e.g. IDT, Sigma, Eurofins):

```tsv
Primer_Name	Sequence_5_to_3	Notes
Uni-Seq-F1	GTGACACATGCCAAAGAA	Universal Sanger Fwd (82.0% cov, 55.26C)
Uni-Seq-R1	AAGAGATGAGGCGACAGT	Universal Sanger Rev (78.0% cov, 57.35C)
Spec_F_Lineage_1	ACAGAAACTGTGGACACG	Strain-Specific (Lineage_1, 2.4%)
Spec_F_Lineage_2	ACAGAAACTGTGGACACA	Strain-Specific (Lineage_2, 1.2%)
```

---
*Report produced by AI4VPrimer Amplicon-Sanger Pipeline.*