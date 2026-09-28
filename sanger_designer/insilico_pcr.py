"""
insilico_pcr.py
Simulates PCR amplification on an aligned sequence dataset using user-provided
or auto-designed PCR primers, validates PCR sensitivity, and extracts the amplicon slice.
"""

from typing import List, Dict, Tuple, Optional
import re
from .msa_utils import reverse_complement, iupac_matches, IUPAC_TABLE, slice_alignment


def iupac_to_regex(iupac_seq: str) -> str:
    """Convert an IUPAC nucleotide sequence into a regex pattern."""
    pattern = []
    for char in iupac_seq.upper():
        bases = IUPAC_TABLE.get(char, {char})
        if len(bases) == 1:
            pattern.append(list(bases)[0])
        else:
            pattern.append("[" + "".join(sorted(bases)) + "]")
    return "".join(pattern)


class InSilicoPCR:
    """
    Simulates in-silico PCR on an alignment to locate primer binding sites,
    calculate PCR amplification sensitivity, and extract the resulting amplicon MSA.
    """

    def __init__(
        self,
        fwd_primer: str,
        rev_primer: str,
        max_mismatches: int = 1,
        require_exact_3prime: int = 3
    ):
        """
        :param fwd_primer: Forward PCR primer sequence (5' to 3')
        :param rev_primer: Reverse PCR primer sequence (5' to 3', user standard orientation)
        :param max_mismatches: Maximum allowable mismatches across the primer
        :param require_exact_3prime: Number of bases at 3' end that must match with 0 mismatches
        """
        self.fwd_primer = fwd_primer.strip().upper()
        self.rev_primer = rev_primer.strip().upper()
        # The reverse primer binds to the top strand in reverse complement orientation
        self.rev_primer_rc = reverse_complement(self.rev_primer)
        self.max_mismatches = max_mismatches
        self.require_exact_3prime = require_exact_3prime

    def _match_primer_to_seq(self, primer: str, template: str, is_reverse: bool = False) -> List[Tuple[int, int, int]]:
        """
        Scan template (without gaps or across gapped positions) for primer binding sites.
        Returns list of (start_idx, end_idx, mismatch_count).
        """
        p_len = len(primer)
        hits = []
        t_len = len(template)

        for i in range(t_len - p_len + 1):
            sub = template[i : i + p_len]
            if '-' in sub or '.' in sub:
                continue
            
            mismatches = 0
            match_possible = True
            
            # Check 3' terminal anchor first
            # For forward primer, 3' is at the end of the window (sub[-3:])
            # For reverse primer annealing, the 3' end on the reverse primer corresponds to sub[:3] of rev_primer_rc
            if not is_reverse:
                # Forward: 3' end is at the right
                for k in range(p_len - self.require_exact_3prime, p_len):
                    if not iupac_matches(primer[k], sub[k]):
                        match_possible = False
                        break
            else:
                # Reverse primer RC: 3' of the original primer is 5' of the RC
                for k in range(self.require_exact_3prime):
                    if not iupac_matches(primer[k], sub[k]):
                        match_possible = False
                        break

            if not match_possible:
                continue

            for j in range(p_len):
                if not iupac_matches(primer[j], sub[j]):
                    mismatches += 1
                    if mismatches > self.max_mismatches:
                        match_possible = False
                        break
            
            if match_possible:
                hits.append((i, i + p_len, mismatches))
                
        return hits

    def run_pcr(
        self,
        headers: List[str],
        sequences: List[str]
    ) -> Dict:
        """
        Run in-silico PCR across all sequences in the alignment.
        Returns a dictionary containing:
          - 'amplified_count': int
          - 'total_count': int
          - 'sensitivity_pct': float
          - 'start_col': consensus start column in MSA
          - 'end_col': consensus end column in MSA
          - 'amplicon_len': length of consensus amplicon
          - 'amplicon_sequences': sliced MSA list
          - 'per_seq_results': detailed per-sequence binding metrics
        """
        num_seqs = len(sequences)
        amplified = 0
        fwd_starts = []
        rev_ends = []
        per_seq = []

        for h, s in zip(headers, sequences):
            f_hits = self._match_primer_to_seq(self.fwd_primer, s, is_reverse=False)
            r_hits = self._match_primer_to_seq(self.rev_primer_rc, s, is_reverse=True)

            valid_pair = None
            if f_hits and r_hits:
                # Find best valid amplicon (fwd before rev)
                for f_start, f_end, f_mis in f_hits:
                    for r_start, r_end, r_mis in r_hits:
                        if f_start < r_end and (r_end - f_start) >= (len(self.fwd_primer) + len(self.rev_primer)):
                            valid_pair = (f_start, r_end, f_mis, r_mis)
                            break
                    if valid_pair:
                        break

            if valid_pair:
                amplified += 1
                fwd_starts.append(valid_pair[0])
                rev_ends.append(valid_pair[1])
                per_seq.append({
                    "id": h,
                    "amplified": True,
                    "start": valid_pair[0],
                    "end": valid_pair[1],
                    "length": valid_pair[1] - valid_pair[0],
                    "fwd_mismatches": valid_pair[2],
                    "rev_mismatches": valid_pair[3]
                })
            else:
                per_seq.append({
                    "id": h,
                    "amplified": False,
                    "fwd_hits": len(f_hits),
                    "rev_hits": len(r_hits)
                })

        sensitivity = (amplified / num_seqs) * 100.0 if num_seqs > 0 else 0.0

        if amplified == 0:
            raise ValueError(
                f"In-Silico PCR failed: No sequence was amplified. "
                f"Verify that Forward ({self.fwd_primer}) and Reverse ({self.rev_primer}) "
                f"match the target gene."
            )

        # Use median/most frequent boundary for consensus amplicon coordinates
        consensus_start = int(sorted(fwd_starts)[len(fwd_starts) // 2])
        consensus_end = int(sorted(rev_ends)[len(rev_ends) // 2])
        amplicon_msa = slice_alignment(sequences, consensus_start, consensus_end)

        return {
            "amplified_count": amplified,
            "total_count": num_seqs,
            "sensitivity_pct": round(sensitivity, 2),
            "start_col": consensus_start,
            "end_col": consensus_end,
            "amplicon_len": consensus_end - consensus_start,
            "amplicon_sequences": amplicon_msa,
            "per_seq_results": per_seq
        }
