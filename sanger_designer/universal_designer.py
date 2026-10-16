"""
universal_designer.py
Designs universal Sanger sequencing primers with:
1. Strict 0-degeneracy (Tier 1) targeting invariant sites (Tm 55-60°C, no >=4 repeats).
2. Dynamic fallback to 2-fold degeneracy (Tier 2) if Tier 1 coverage < 80%,
   enforcing the strict 3'-End Protection Rule (no wobble in 3' last 5 nt)
   and attaching a mandatory '2x Concentration' lab warning.
"""

from typing import List, Dict, Optional, Tuple
import re
from .msa_utils import (
    get_column_frequencies,
    reverse_complement,
    iupac_matches,
    IUPAC_TWO_FOLD,
    IUPAC_TABLE
)

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


def calculate_tm_simple(seq_str: str) -> float:
    """Basic nearest-neighbor / GC-based Tm approximation if primer3 is unavailable."""
    seq = seq_str.upper()
    g = seq.count('G') + seq.count('S') + 0.5 * (seq.count('R') + seq.count('Y') + seq.count('K') + seq.count('M'))
    c = seq.count('C') + seq.count('S') + 0.5 * (seq.count('R') + seq.count('Y') + seq.count('K') + seq.count('M'))
    gc_count = g + c
    length = len(seq)
    if length == 0:
        return 0.0
    # Standard formula for 14-70 nt: 64.9 + 41 * (yG + zC - 16.4) / length
    return round(64.9 + 41.0 * (gc_count - 16.4) / length, 2)


def expand_iupac_variants(seq_str: str) -> List[str]:
    """Expands an IUPAC degenerate sequence into all unambiguous ACGT variants."""
    variants = [""]
    for char in seq_str.upper():
        allowed = IUPAC_TABLE.get(char, {char})
        variants = [v + b for v in variants for b in allowed]
    return variants


def check_homopolymers(seq_str: str, max_repeat: int = 3) -> bool:
    """Returns True if any nucleotide is repeated more than max_repeat consecutive times."""
    pattern = r"(A{" + str(max_repeat + 1) + r",}|C{" + str(max_repeat + 1) + r",}|G{" + str(max_repeat + 1) + r",}|T{" + str(max_repeat + 1) + r",})"
    return bool(re.search(pattern, seq_str.upper()))


def evaluate_biophysics(seq_str: str, min_tm: float = 55.0, max_tm: float = 60.0) -> Dict:
    """Calculates Tm, GC%, homopolymers, and secondary structure stability."""
    seq = seq_str.upper()
    length = len(seq)
    gc_count = sum(1 for b in seq if b in ('G', 'C', 'S'))
    gc_pct = round((gc_count / length) * 100.0, 2)

    has_homopolymer = check_homopolymers(seq, max_repeat=3)

    tm = 0.0
    hairpin_dg = 0.0
    homodimer_dg = 0.0
    passed_tm = False

    if HAS_PRIMER3:
        variants = expand_iupac_variants(seq)
        tms = []
        hairpin_dgs = []
        homodimer_dgs = []
        for var in variants:
            tms.append(primer3.calc_tm(var, mv_conc=50, dv_conc=1.5, dna_conc=200))
            try:
                hp = primer3.calc_hairpin(var)
                hairpin_dgs.append(hp.dg / 1000.0)
            except Exception:
                hairpin_dgs.append(0.0)
            try:
                hd = primer3.calc_homodimer(var)
                homodimer_dgs.append(hd.dg / 1000.0)
            except Exception:
                homodimer_dgs.append(0.0)
        tm = round(sum(tms) / len(tms), 2)
        hairpin_dg = round(min(hairpin_dgs), 2)
        homodimer_dg = round(min(homodimer_dgs), 2)
        passed_tm = all(min_tm <= t <= max_tm for t in tms)
    else:
        tm = calculate_tm_simple(seq)
        passed_tm = (min_tm <= tm <= max_tm)

    passed_gc = (40.0 <= gc_pct <= 60.0)
    passed_repeat = not has_homopolymer
    passed_secondary = (hairpin_dg > -3.0 and homodimer_dg > -4.5)

    is_valid = passed_tm and passed_gc and passed_repeat and passed_secondary

    return {
        "sequence": seq,
        "length": length,
        "tm": tm,
        "gc_pct": gc_pct,
        "has_homopolymer": has_homopolymer,
        "hairpin_dg": hairpin_dg,
        "homodimer_dg": homodimer_dg,
        "is_valid": is_valid
    }


class UniversalSangerDesigner:
    """
    Designs Universal Sanger sequencing primers with a 2-Tier discovery architecture.
    """

    def __init__(
        self,
        min_coverage_pct: float = 80.0,
        min_len: int = 18,
        max_len: int = 24,
        min_tm: float = 55.0,
        max_tm: float = 60.0,
        ideal_tm: float = 57.5,
        dye_blob_buffer: int = 40
    ):
        self.min_coverage_pct = min_coverage_pct
        self.min_len = min_len
        self.max_len = max_len
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.ideal_tm = ideal_tm
        self.dye_blob_buffer = dye_blob_buffer

    def _calculate_sequence_coverage(self, primer_seq: str, sequences: List[str], start_col: int, end_col: int) -> float:
        """Percentage of sequences in alignment that match the primer sequence with 0 mismatches."""
        num_seqs = len(sequences)
        if num_seqs == 0:
            return 0.0
        matches = 0
        p_len = len(primer_seq)

        # High-performance C-speed fast path for invariant (Tier 1) non-degenerate primers
        is_pure_dna = all(c in "ACGTacgt" for c in primer_seq)
        if is_pure_dna:
            for s in sequences:
                if s[start_col:end_col] == primer_seq:
                    matches += 1
            return round((matches / num_seqs) * 100.0, 2)

        # High-performance set lookup for 2-fold degenerate primers
        variants_set = set(expand_iupac_variants(primer_seq))
        for s in sequences:
            sub = s[start_col:end_col]
            if sub in variants_set:
                matches += 1
            elif '-' in sub or '.' in sub or len(sub) != p_len:
                continue
            elif any(c not in "ACGT" for c in sub):
                if all(iupac_matches(primer_seq[k], sub[k]) for k in range(p_len)):
                    matches += 1

        return round((matches / num_seqs) * 100.0, 2)

    def design_universal_primer(
        self,
        amplicon_msa: List[str],
        direction: str = "FORWARD",
        target_subregion_start: Optional[int] = None,
        target_subregion_end: Optional[int] = None,
        exclude_range: Optional[Tuple[int, int]] = None,
        col_freqs: Optional[List[Dict[str, float]]] = None,
        search_window: Optional[Tuple[int, int]] = None
    ) -> Dict:
        """
        Designs the optimal universal Sanger sequencing primer.
        :param amplicon_msa: Aligned sequences of the amplicon.
        :param direction: "FORWARD" or "REVERSE".
        :param target_subregion_start: Start index of the target subregion (e.g. variable region) relative to amplicon.
        :param target_subregion_end: End index of target subregion.
        :param exclude_range: Optional (start, end) column range to strictly prevent forward/reverse overlap.
        :param col_freqs: Optional precalculated column frequencies.
        :param search_window: Optional explicit (min_pos, max_pos) tuple.
        :return: Dict containing selected candidate, tier level, concentration instructions, and runner-ups.
        """
        amp_len = len(amplicon_msa[0])
        num_seqs = len(amplicon_msa)
        if col_freqs is None:
            col_freqs = get_column_frequencies(amplicon_msa)

        # Precompute per-column consensus metrics
        min_cov_ratio = self.min_coverage_pct / 100.0

        top1_bases = []
        top1_freqs = []
        top2_bases = []
        top2_iupacs = []
        top2_freqs = []

        for c in range(amp_len):
            freq = col_freqs[c]
            if not freq:
                top1_bases.append('N')
                top1_freqs.append(0.0)
                top2_bases.append(None)
                top2_iupacs.append(None)
                top2_freqs.append(0.0)
                continue

            sorted_b = sorted(freq.items(), key=lambda x: x[1], reverse=True)
            top1_bases.append(sorted_b[0][0])
            top1_freqs.append(sorted_b[0][1])

            if len(sorted_b) >= 2:
                top2_freqs.append(sorted_b[0][1] + sorted_b[1][1])
                pair = frozenset([sorted_b[0][0], sorted_b[1][0]])
                top2_bases.append((sorted_b[0][0], sorted_b[1][0]))
                top2_iupacs.append(IUPAC_TWO_FOLD.get(pair))
            else:
                top2_freqs.append(sorted_b[0][1])
                top2_bases.append(None)
                top2_iupacs.append(None)

        # Determine search window boundaries
        if search_window is not None:
            min_pos, max_pos = search_window
            min_pos = max(0, min(amp_len - self.min_len, min_pos))
            max_pos = max(min_pos + self.min_len, min(amp_len, max_pos))
        elif direction.upper() == "FORWARD":
            if target_subregion_start is not None:
                max_pos = max(0, target_subregion_start - self.dye_blob_buffer)
                min_pos = max(0, max_pos - 250)
            else:
                min_pos = 0 if amp_len <= 500 else 35
                max_pos = amp_len // 2 if amp_len <= 500 else min(amp_len // 2, 350)
        else: # REVERSE primer: must sit DOWNSTREAM of target_subregion or at 3' terminal section
            if target_subregion_end is not None:
                min_pos = min(amp_len, target_subregion_end + self.dye_blob_buffer)
                max_pos = min(amp_len, min_pos + 250)
            else:
                min_pos = amp_len // 2 if amp_len <= 500 else max(amp_len // 2, amp_len - 350)
                max_pos = amp_len if amp_len <= 500 else max(min_pos + 1, amp_len - 35)

        # ----------------------------------------------------
        # TIER 1: STRICT INVARIANT (0 Degeneracy)
        # ----------------------------------------------------
        tier1_candidates = []

        for p_len in range(self.min_len, self.max_len + 1):
            for start_idx in range(min_pos, max_pos - p_len + 1):
                end_idx = start_idx + p_len

                # Enforce exclusion range (e.g. no overlap with opposing primer)
                if exclude_range:
                    ex_start, ex_end = exclude_range
                    if not (end_idx <= ex_start or start_idx >= ex_end):
                        continue

                # CONSENSUS FILTER FOR TIER 1:
                # If ANY position has top base frequency < min_cov_ratio (80%),
                # the unambiguous primer can never achieve 80% coverage. Skip immediately!
                if any(top1_freqs[c] < min_cov_ratio for c in range(start_idx, end_idx)):
                    continue

                # Build consensus oligo from top base at each column
                candidate_seq = "".join(top1_bases[c] for c in range(start_idx, end_idx))
                if 'N' in candidate_seq or '-' in candidate_seq:
                    continue

                final_seq = candidate_seq if direction.upper() == "FORWARD" else reverse_complement(candidate_seq)

                # Fast biophysical pre-filter:
                # 1. Homopolymers check
                if check_homopolymers(final_seq, max_repeat=3):
                    continue

                # 2. Fast GC% check (40.0 - 60.0%)
                gc_count = sum(1 for b in final_seq if b in ('G', 'C', 'S'))
                gc_pct = (gc_count / p_len) * 100.0
                if not (40.0 <= gc_pct <= 60.0):
                    continue

                # 3. Primer3 biophysical validation (Tm, Hairpin, Homodimer)
                bio = evaluate_biophysics(final_seq, self.min_tm, self.max_tm)
                if not bio["is_valid"]:
                    continue

                # 4. Alignment sequence coverage calculation
                coverage = self._calculate_sequence_coverage(
                    candidate_seq, amplicon_msa, start_idx, end_idx
                )

                if coverage >= self.min_coverage_pct:
                    tier1_candidates.append({
                        "tier": 1,
                        "sequence": final_seq,
                        "order_sequence": final_seq,
                        "direction": direction.upper(),
                        "start_col": start_idx,
                        "end_col": end_idx,
                        "length": p_len,
                        "tm": bio["tm"],
                        "gc_pct": bio["gc_pct"],
                        "hairpin_dg": bio["hairpin_dg"],
                        "homodimer_dg": bio["homodimer_dg"],
                        "coverage_pct": coverage,
                        "degeneracy_fold": 1,
                        "requires_2x_concentration": False,
                        "concentration_recommendation": "Standard 1x (e.g. 3.2 pmol per 10 uL reaction)",
                        "notes": "Optimal 0-degeneracy invariant primer."
                    })

        if tier1_candidates:
            tier1_candidates.sort(key=lambda x: (-x["coverage_pct"], abs(x["tm"] - self.ideal_tm)))
            return {
                "status": "SUCCESS",
                "tier_selected": 1,
                "selected_primer": tier1_candidates[0],
                "all_candidates": tier1_candidates[:5]
            }

        # ----------------------------------------------------
        # TIER 2: FALLBACK (2-Fold Degeneracy)
        # Triggered when Tier 1 coverage is below threshold (<80%)
        # ----------------------------------------------------
        tier2_candidates = []

        for p_len in range(self.min_len, self.max_len + 1):
            for start_idx in range(min_pos, max_pos - p_len + 1):
                end_idx = start_idx + p_len

                # Enforce exclusion range in Tier 2 as well
                if exclude_range:
                    ex_start, ex_end = exclude_range
                    if not (end_idx <= ex_start or start_idx >= ex_end):
                        continue

                # MASTER CONSENSUS PRE-FILTER FOR TIER 2:
                # Count positions with low consensus (< min_cov_ratio)
                low_consensus_positions = [
                    r for r, c in enumerate(range(start_idx, end_idx))
                    if top1_freqs[c] < min_cov_ratio
                ]

                # PRINCIPLE RULE: If >= 2 positions have low consensus,
                # DISCARD WINDOW OF THIS 20nt IMMEDIATELY!
                # (Because 1 degenerate base can only rescue 1 position)
                if len(low_consensus_positions) >= 2:
                    continue

                # Identify the single wobble target position
                if len(low_consensus_positions) == 1:
                    wobble_rel_pos = low_consensus_positions[0]
                else:
                    # len(low_consensus_positions) == 0:
                    # Target the single position with lowest top1 frequency
                    # that satisfies the 3'-end protection rule
                    valid_wobble_positions = [
                        r for r in range(p_len)
                        if (r < (p_len - 5) if direction.upper() == "FORWARD" else r >= 5)
                    ]
                    if not valid_wobble_positions:
                        continue
                    wobble_rel_pos = min(valid_wobble_positions, key=lambda r: top1_freqs[start_idx + r])

                # STRICT 3'-END PROTECTION RULE:
                # No wobble allowed within the last 5 nucleotides of the 3' end!
                if direction.upper() == "FORWARD":
                    if wobble_rel_pos >= (p_len - 5):
                        continue
                else: # REVERSE primer (3' end is at relative 0..4 of top strand)
                    if wobble_rel_pos < 5:
                        continue

                wobble_col = start_idx + wobble_rel_pos

                # Check top-2 base combined frequency >= min_cov_ratio
                if top2_freqs[wobble_col] < min_cov_ratio:
                    continue

                iupac_char = top2_iupacs[wobble_col]
                if not iupac_char:
                    continue

                # Construct degenerate primer sequence
                candidate_chars = []
                for r, c in enumerate(range(start_idx, end_idx)):
                    if r == wobble_rel_pos:
                        candidate_chars.append(iupac_char)
                    else:
                        candidate_chars.append(top1_bases[c])

                if 'N' in candidate_chars or '-' in candidate_chars:
                    continue

                top_strand_seq = "".join(candidate_chars)
                final_seq = top_strand_seq if direction.upper() == "FORWARD" else reverse_complement(top_strand_seq)

                # Fast biophysical pre-filter
                if check_homopolymers(final_seq, max_repeat=3):
                    continue

                gc_count = sum(1 for b in final_seq if b in ('G', 'C', 'S'))
                gc_pct = (gc_count / p_len) * 100.0
                if not (40.0 <= gc_pct <= 60.0):
                    continue

                # Primer3 biophysical validation
                bio = evaluate_biophysics(final_seq, self.min_tm, self.max_tm)
                if not bio["is_valid"]:
                    continue

                # Coverage calculation across alignment
                coverage = self._calculate_sequence_coverage(
                    top_strand_seq, amplicon_msa, start_idx, end_idx
                )

                if coverage >= self.min_coverage_pct:
                    tier2_candidates.append({
                        "tier": 2,
                        "sequence": final_seq,
                        "order_sequence": final_seq,
                        "direction": direction.upper(),
                        "start_col": start_idx,
                        "end_col": end_idx,
                        "length": p_len,
                        "tm": bio["tm"],
                        "gc_pct": bio["gc_pct"],
                        "hairpin_dg": bio["hairpin_dg"],
                        "homodimer_dg": bio["homodimer_dg"],
                        "coverage_pct": coverage,
                        "degeneracy_fold": 2,
                        "wobble_position_from_5prime": wobble_rel_pos + 1 if direction.upper() == "FORWARD" else p_len - wobble_rel_pos,
                        "requires_2x_concentration": True,
                        "concentration_recommendation": (
                            "⚠️ Double primer concentration: Use 6.4 pmol per 10 uL reaction "
                            "(2x standard) to compensate for 50% signal dilution from the 2-fold degenerate oligo."
                        ),
                        "notes": "Tier 2 Fallback: 2-fold degenerate base safely placed outside the 3' anchor zone."
                    })

        if tier2_candidates:
            tier2_candidates.sort(key=lambda x: (-x["coverage_pct"], abs(x["tm"] - self.ideal_tm)))
            return {
                "status": "SUCCESS",
                "tier_selected": 2,
                "selected_primer": tier2_candidates[0],
                "all_candidates": tier2_candidates[:5]
            }

        return {
            "status": "FAILED",
            "tier_selected": None,
            "selected_primer": None,
            "message": (
                f"No universal primer met the {self.min_coverage_pct}% coverage threshold, "
                "even with 2-fold degeneracy fallback. Please use Mode 2 (Strain-Specific Panel) "
                "to sequence individual variants in separate tubes."
            )
        }
