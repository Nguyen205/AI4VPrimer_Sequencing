"""
test_pipeline.py
Unit tests and validation suite for the Amplicon Sanger Sequencing Primer Design pipeline.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from sanger_designer.msa_utils import (
    reverse_complement,
    iupac_matches,
    load_alignment,
    IUPAC_TWO_FOLD
)
from sanger_designer.insilico_pcr import InSilicoPCR
from sanger_designer.universal_designer import UniversalSangerDesigner, evaluate_biophysics
from sanger_designer.strain_specific_designer import StrainSpecificSangerDesigner
from sanger_designer.tiling import SangerTiler
from sanger_designer.pipeline import SangerAmpliconPipeline


class TestSangerPrimerDesigner(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        os.makedirs(self.test_dir, exist_ok=True)
        self.fasta_path = os.path.join(self.test_dir, "synthetic_test.fasta")

        # Construct a synthetic aligned test dataset of 10 viral strains (length 600 bp)
        # - Position 0-25: PCR Forward binding site (invariant: "ATGCGATCGATCGATCGATCGA")
        # - Position 50-70: Conserved region (high Tm ~57C, 100% invariant)
        # - Position 200-220: Variable region with 2 distinct haplotypes (Group A: 6 seqs, Group B: 4 seqs)
        # - Position 400-420: Diverse region where an A/G transition occurs at pos 405 (eligible for 2-fold degenerate)
        # - Position 575-600: PCR Reverse binding site
        fwd_pcr = "ATGCGATCGATCGATCGATCGA" # 22 nt
        rev_pcr = "CGTAACGTAGCTAGCTAGCTA"  # 21 nt
        rev_pcr_rc = reverse_complement(rev_pcr)

        seqs = []
        mixed_spacer = "AGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGATCGTAGCTA"
        for i in range(10):
            prefix = fwd_pcr + (mixed_spacer[:28]) # pos 0 to 50
            # Conserved Sanger Fwd site (pos 50 to 72, 22nt, Tm ~57C)
            sanger_fwd_site = "GGATCCGAATGCTGTAGCACCA" 
            spacer1 = (mixed_spacer * 3)[:128]
            # Variable / Discriminatory site at pos 200:
            # 6 seqs have ...ACCG, 4 seqs have ...ACCT (3' mismatch between them!)
            var_site = "TGCAAGTTCCGAACGAGACC" + ("G" if i < 6 else "T")
            spacer2 = (mixed_spacer * 4)[:179]
            # Region with 2-fold wobble (A for 5 seqs, G for 5 seqs at pos 4)
            wobble_site = "CTAG" + ("A" if i < 5 else "G") + "TGCGAATGCTGACCG"
            spacer3 = (mixed_spacer * 4)[:154]
            full_seq = prefix + sanger_fwd_site + spacer1 + var_site + spacer2 + wobble_site + spacer3 + rev_pcr_rc
            seqs.append((f">strain_{i+1}", full_seq[:650]))

        with open(self.fasta_path, "w") as f:
            for header, s in seqs:
                f.write(f"{header}\n{s}\n")

        self.fwd_pcr = fwd_pcr
        self.rev_pcr = rev_pcr

    def test_reverse_complement_and_iupac(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AYGR"), "YCRT")
        self.assertTrue(iupac_matches("R", "A"))
        self.assertTrue(iupac_matches("R", "G"))
        self.assertFalse(iupac_matches("R", "C"))
        self.assertTrue(iupac_matches("Y", "C"))

    def test_insilico_pcr(self):
        headers, sequences = load_alignment(self.fasta_path)
        pcr = InSilicoPCR(self.fwd_pcr, self.rev_pcr)
        res = pcr.run_pcr(headers, sequences)
        self.assertEqual(res["amplified_count"], 10)
        self.assertEqual(res["sensitivity_pct"], 100.0)
        self.assertGreater(res["amplicon_len"], 500)

    def test_universal_sanger_designer_tier1(self):
        headers, sequences = load_alignment(self.fasta_path)
        designer = UniversalSangerDesigner(min_coverage_pct=80.0, min_tm=50.0, max_tm=65.0)
        res = designer.design_universal_primer(sequences, direction="FORWARD")
        self.assertEqual(res["status"], "SUCCESS")
        self.assertEqual(res["tier_selected"], 1)
        self.assertGreaterEqual(res["selected_primer"]["coverage_pct"], 80.0)
        self.assertFalse(res["selected_primer"]["requires_2x_concentration"])

    def test_strain_specific_panel(self):
        headers, sequences = load_alignment(self.fasta_path)
        designer = StrainSpecificSangerDesigner(min_tm=50.0, max_tm=65.0)
        # Ask to design forward primers upstream of variable region (pos 200)
        panel_res = designer.design_panel(sequences, headers, direction="FORWARD", target_subregion_start=205)
        self.assertEqual(panel_res["status"], "SUCCESS")
        self.assertGreaterEqual(panel_res["panel_size"], 1)
        self.assertGreaterEqual(panel_res["total_coverage_pct"], 60.0)

    def test_tiling_thresholds(self):
        tiler = SangerTiler(tiling_threshold=1200)
        headers, sequences = load_alignment(self.fasta_path)
        designer = UniversalSangerDesigner(min_tm=50.0, max_tm=65.0)
        
        # Test 1: 650 bp amplicon -> No tiling
        res1 = tiler.evaluate_and_tile(sequences, designer)
        self.assertFalse(res1["tiling_required"])

        # Test 2: 1500 bp synthetic amplicon -> Tiling activated
        long_seqs = [s * 3 for s in sequences]
        res2 = tiler.evaluate_and_tile(long_seqs, designer)
        self.assertTrue(res2["tiling_required"])
        self.assertGreater(len(res2["tiled_primers"]), 1)

    def test_full_pipeline_run(self):
        report_out = os.path.join(self.test_dir, "test_report.md")
        pipeline = SangerAmpliconPipeline(
            fasta_path=self.fasta_path,
            fwd_pcr_primer=self.fwd_pcr,
            rev_pcr_primer=self.rev_pcr,
            min_tm=50.0,
            max_tm=65.0,
            output_report_path=report_out
        )
        res = pipeline.run()
        self.assertIn("report_content", res)
        self.assertTrue(os.path.exists(report_out))
        with open(report_out) as f:
            content = f.read()
        self.assertIn("PCR Amplicon Summary", content)
        self.assertIn("Mode 1: Universal Sanger Sequencing Primers", content)
        self.assertIn("Mode 2: Strain-Specific Panel", content)


if __name__ == "__main__":
    unittest.main()
