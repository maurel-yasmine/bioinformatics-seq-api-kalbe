from typing import Dict, List
from app.services.translator import STANDARD_GENETIC_CODE


def annotate_variants(reference_seq: str, variants: List[Dict]) -> List[Dict]:
    annotated: List[Dict] = []

    for var in variants:
        v_type = var["type"]
        pos = var["position"]
        ref_val = var["ref"]
        alt_val = var["alt"]
        effect = "undetermined"

        if v_type == "SNV":
            codon_index = (pos - 1) // 3
            codon_start = codon_index * 3
            if codon_start + 3 <= len(reference_seq):
                ref_codon = reference_seq[codon_start : codon_start + 3]
                offset = (pos - 1) % 3
                alt_codon = ref_codon[:offset] + alt_val + ref_codon[offset + 1 :]

                if "N" not in ref_codon and "N" not in alt_codon:
                    ref_aa = STANDARD_GENETIC_CODE[ref_codon]
                    alt_aa = STANDARD_GENETIC_CODE[alt_codon]
                    if alt_aa == "*" and ref_aa != "*":
                        effect = "nonsense"
                    elif ref_aa == alt_aa:
                        effect = "synonymous"
                    else:
                        effect = "missense"

        elif v_type == "deletion":
            effect = "in-frame deletion" if len(ref_val) % 3 == 0 else "frameshift"

        elif v_type == "insertion":
            effect = "in-frame insertion" if len(alt_val) % 3 == 0 else "frameshift"

        annotated.append(
            {
                "position": pos,
                "ref": ref_val,
                "alt": alt_val,
                "type": v_type,
                "effect": effect,
            }
        )

    return annotated