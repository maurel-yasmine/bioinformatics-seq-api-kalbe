from typing import Dict, List, Tuple

STANDARD_GENETIC_CODE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_sequence(sequence: str, frame: int = 1) -> Tuple[str, int, bool]:
    if frame not in (1, 2, 3):
        raise ValueError("Frame must be 1, 2, or 3.")

    coding_seq = sequence[frame - 1 :]
    amino_acids: List[str] = []
    has_stop_codon = False

    for i in range(0, len(coding_seq) - 2, 3):
        codon = coding_seq[i : i + 3]
        if "N" in codon:
            amino_acids.append("X")
            continue

        aa = STANDARD_GENETIC_CODE[codon]
        if aa == "*":
            has_stop_codon = True
            break
        amino_acids.append(aa)

    protein = "".join(amino_acids)
    return protein, len(protein), has_stop_codon