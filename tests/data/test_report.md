# Amplicon & Sanger Sequencing Primer Design Report

Generated automatically by **AI4VPrimer Amplicon-Sanger Engine**.

---

## 1. PCR Amplicon Summary
- **Mode:** External PCR Primers Provided (Path A)
- **Template FASTA:** `/Users/aijiazhou/Desktop/Amplicon_Sanger_Primer_Suite/tests/data/synthetic_test.fasta` (10 sequences)
- **Forward PCR Primer:** `ATGCGATCGATCGATCGATCGA`
- **Reverse PCR Primer:** `CGTAACGTAGCTAGCTAGCTA`
- **In-Silico PCR Sensitivity:** **100.0%** (10/10 sequences amplified)
- **Amplicon Span:** MSA columns `0` to `595`
- **Amplicon Length:** **595 bp**

---

## 2. Amplicon Sizing & Tiling Strategy
- **Amplicon Length:** 595 bp
- **Tiling Required:** `NO (<= 1200 bp)`
- **Strategy:** Single / Terminal Forward & Reverse Reactions (No internal walking needed)
- **Details:** Amplicon length (595 bp) is within the <= 1200 bp standard Sanger limit. Opposing Forward and Reverse sequencing reactions overlap in the center to cover the full length.

---

## 3. Mode 1: Universal Sanger Sequencing Primers (Single-Strain Testing)
Use these primers first. If your sample contains a single pure strain, these will yield clean, sharp chromatograms.

### Universal Primers Summary Table
| Primer ID | Role | Direction | Position (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Population Coverage | Tier |
|---|---|---|---|---|---|---|---|---|---|
| `Uni-Seq-F1` | Forward Terminal | FORWARD | **Cols 197 – 215** | 18 nt | `ATCTGCAAGTTCCGAACG` | 57.52 | 50.0% | **100.0%** (10/10 seqs) | Tier 1 (Invariant) |
| `Uni-Seq-R1` | Reverse Terminal | REVERSE | **Cols 413 – 431** | 18 nt | `ATCGATCAGCTCGGTCAG` | 58.36 | 55.56% | **100.0%** (10/10 seqs) | Tier 1 (Invariant) |

### Forward Sequencing Primer (`Uni-Seq-F1`)
- **Tier:** 🟢 Tier 1 (Strict Invariant, 0 Degeneracy)
- **Binding Position:** **Alignment Columns `197 – 215`**
- **Sequence (5' → 3'):** `ATCTGCAAGTTCCGAACG`
- **Length:** 18 nt | **Tm:** 57.52 °C | **GC:** 50.0%
- **Population Coverage:** **100.0%** (10 / 10 strains in `synthetic_test.fasta`)
- **Secondary Structures:** Hairpin $\Delta G$: 0.05 kcal/mol | Dimer $\Delta G$: -3.83 kcal/mol
- **Concentration Protocol:** Standard 1x (e.g. 3.2 pmol per 10 uL reaction)

### Reverse Sequencing Primer (`Uni-Seq-R1`)
- **Tier:** 🟢 Tier 1 (Strict Invariant, 0 Degeneracy)
- **Binding Position:** **Alignment Columns `413 – 431`**
- **Sequence (5' → 3'):** `ATCGATCAGCTCGGTCAG`
- **Length:** 18 nt | **Tm:** 58.36 °C | **GC:** 55.56%
- **Population Coverage:** **100.0%** (10 / 10 strains in `synthetic_test.fasta`)
- **Secondary Structures:** Hairpin $\Delta G$: -0.88 kcal/mol | Dimer $\Delta G$: -3.42 kcal/mol
- **Concentration Protocol:** Standard 1x (e.g. 3.2 pmol per 10 uL reaction)

---

## 4. Mode 2: Strain-Specific Panel (Mixture Deconvolution)
If your clinical/field sample is a **mixture of strains**, universal primers will produce double peaks. 
Use the minimal panel below. Each primer selectively amplifies its specific clade via **3'-terminal extension termination**.

### How Lineages are Determined
1. **Target Region Partitioning:** The alignment is clustered by shared nucleotide motifs across the variable target domain, grouping identical alleles into biological lineages.
2. **3'-Terminal Selectivity:** Primers are placed so the very last 3' nucleotide matches only the intended lineage with 0 mismatches, but mismatches all other lineages, arresting non-specific Taq polymerase extension.

- **Clusters/Variants Discovered:** 3
- **Panel Tubes Required:** 3
- **Total Cumulative Coverage:** **240.0%** across all 10 sequences

### Lineage Breakdown
| Lineage ID | Sequence Count | Population Share | Assigned Test Tube |
|---|---|---|---|
| **Lineage_1** | 10 / 10 seqs | **100.0%** | **Tube_1** (`Spec_F_Lineage_1`) |
| **Lineage_2** | 4 / 10 seqs | **40.0%** | **Tube_2** (`Spec_F_Lineage_2`) |
| **Lineage_3** | 10 / 10 seqs | **100.0%** | **Tube_3** (`Spec_F_Lineage_3`) |

### Minimal Panel Tubes
| Tube ID | Primer ID | Target Variant | Position (bp) | Direction | Sequence (5' → 3') | Length | Tm (°C) | Population Share | Reaction Protocol |
|---|---|---|---|---|---|---|---|---|---|
| **Tube_1** | `Spec_F_Lineage_1` | Lineage_1 | **Cols 50 – 68** | FORWARD | `GGATCCGAATGCTGTAGC` | 18 nt | 57.43 | **100.0%** (10 seqs) | Test in Tube 1. Selective 3' extension ensures clean trace for Lineage_1 even in mixtures. |
| **Tube_2** | `Spec_F_Lineage_2` | Lineage_2 | **Cols 203 – 221** | FORWARD | `AAGTTCCGAACGAGACCT` | 18 nt | 57.61 | **40.0%** (4 seqs) | Test in Tube 2. Selective 3' extension ensures clean trace for Lineage_2 even in mixtures. |
| **Tube_3** | `Spec_F_Lineage_3` | Lineage_3 | **Cols 50 – 68** | FORWARD | `GGATCCGAATGCTGTAGC` | 18 nt | 57.43 | **100.0%** (10 seqs) | Test in Tube 3. Selective 3' extension ensures clean trace for Lineage_3 even in mixtures. |

---

## 5. Oligo Synthesis Order Sheet
Copy and paste these sequences directly into your synthesis provider (e.g. IDT, Sigma, Eurofins):

```tsv
Primer_Name	Sequence_5_to_3	Position	Length	Tm_C	Notes
Uni-Seq-F1	ATCTGCAAGTTCCGAACG	Cols 197-215	18	57.52	Universal Sanger Fwd (100.0% cov)
Uni-Seq-R1	ATCGATCAGCTCGGTCAG	Cols 413-431	18	58.36	Universal Sanger Rev (100.0% cov)
Spec_F_Lineage_1	GGATCCGAATGCTGTAGC	Cols 50-68	18	57.43	Strain-Specific (Lineage_1, 100.0%)
Spec_F_Lineage_2	AAGTTCCGAACGAGACCT	Cols 203-221	18	57.61	Strain-Specific (Lineage_2, 40.0%)
Spec_F_Lineage_3	GGATCCGAATGCTGTAGC	Cols 50-68	18	57.43	Strain-Specific (Lineage_3, 100.0%)
```

---
*Report produced by AI4VPrimer Amplicon-Sanger Pipeline.*