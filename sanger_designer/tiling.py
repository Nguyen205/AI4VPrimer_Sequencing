"""
tiling.py
Handles amplicon sizing logic and Sanger primer walking/tiling.
Rule:
- If amplicon <= 1200 bp: Tiling is NOT needed (terminal Forward and Reverse primers suffice).
- If amplicon > 1200 bp: Automatically calculates internal walking steps (~500-650 bp intervals)
  and designs chained sequencing primers to guarantee complete end-to-end coverage.
"""

from typing import List, Dict, Optional
from .universal_designer import UniversalSangerDesigner
from .msa_utils import get_column_frequencies


class SangerTiler:
    """
    Manages amplicon tiling and primer walking.
    """

    def __init__(
        self,
        tiling_threshold: int = 1200,
        walking_step_size: int = 550,
        sanger_read_span: int = 700
    ):
        self.tiling_threshold = tiling_threshold
        self.walking_step_size = walking_step_size
        self.sanger_read_span = sanger_read_span

    def evaluate_and_tile(
        self,
        amplicon_msa: List[str],
        universal_designer: UniversalSangerDesigner,
        col_freqs: Optional[List[Dict[str, float]]] = None
    ) -> Dict:
        """
        Determines if tiling is required and generates tiled sequencing primers if amplicon > 1200 bp.
        """
        amp_len = len(amplicon_msa[0])
        if col_freqs is None:
            col_freqs = get_column_frequencies(amplicon_msa)

        is_tiling_required = (amp_len > self.tiling_threshold)

        if not is_tiling_required:
            # Single-step terminal reactions suffice
            fwd_res = universal_designer.design_universal_primer(
                amplicon_msa, direction="FORWARD", col_freqs=col_freqs
            )
            rev_res = universal_designer.design_universal_primer(
                amplicon_msa, direction="REVERSE", col_freqs=col_freqs
            )
            return {
                "amplicon_length": amp_len,
                "tiling_required": False,
                "strategy": "Single / Terminal Forward & Reverse Reactions (No internal walking needed)",
                "explanation": (
                    f"Amplicon length ({amp_len} bp) is within the <= 1200 bp standard Sanger limit. "
                    "Opposing Forward and Reverse sequencing reactions overlap in the center to cover the full length."
                ),
                "forward_primer": fwd_res.get("selected_primer"),
                "reverse_primer": rev_res.get("selected_primer"),
                "tiled_primers": []
            }

        # Amplicon > 1200 bp: Walking / Tiling algorithm activated
        tiled_primers = []
        step = 0
        current_pos = 0

        while current_pos < (amp_len - 300):
            step += 1
            if step == 1:
                search_win = (35, min(350, amp_len // 2))
            else:
                search_win = (max(0, current_pos - 150), min(amp_len - 35, current_pos + 150))

            sub_design = universal_designer.design_universal_primer(
                amplicon_msa,
                direction="FORWARD",
                search_window=search_win,
                col_freqs=col_freqs
            )

            primer_data = sub_design.get("selected_primer")
            if primer_data:
                primer_data["tile_step"] = step
                primer_data["covered_range"] = f"{primer_data['start_col']} - {min(amp_len, primer_data['start_col'] + self.sanger_read_span)} bp"
                tiled_primers.append(primer_data)
                # Next step
                current_pos = primer_data["start_col"] + self.walking_step_size
            else:
                # If no primer found at this exact window, advance slightly
                current_pos += 150

        # Always ensure a terminal reverse primer is included if found
        rev_res = universal_designer.design_universal_primer(
            amplicon_msa, direction="REVERSE", col_freqs=col_freqs
        )
        if rev_res.get("selected_primer"):
            rev_p = rev_res["selected_primer"]
            rev_p["tile_step"] = len(tiled_primers) + 1
            rev_p["covered_range"] = f"{max(0, rev_p['end_col'] - self.sanger_read_span)} - {rev_p['end_col']} bp (Reverse)"
            tiled_primers.append(rev_p)

        return {
            "amplicon_length": amp_len,
            "tiling_required": True,
            "strategy": f"Multi-Primer Walking ({len(tiled_primers)} primers spaced every ~{self.walking_step_size} bp)",
            "explanation": (
                f"Amplicon length ({amp_len} bp) exceeds the 1200 bp threshold. "
                f"{len(tiled_primers)} chained sequencing primers designed to achieve full-length coverage."
            ),
            "tiled_primers": tiled_primers
        }
