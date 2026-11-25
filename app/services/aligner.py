from typing import Dict, List, Tuple

MATCH_SCORE = 1.0
MISMATCH_SCORE = -1.0
AMBIGUOUS_N_SCORE = 0.0
GAP_PENALTY = -2.0


def _score_pair(ref_base: str, query_base: str) -> float:
    if ref_base == "N" or query_base == "N":
        return AMBIGUOUS_N_SCORE
    if ref_base == query_base:
        return MATCH_SCORE
    return MISMATCH_SCORE


def needleman_wunsch(reference: str, query: str) -> Tuple[float, str, str]:
    m = len(reference)
    n = len(query)

    dp = [[0.0] * (n + 1) for _ in range(m + 1)]
    traceback = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = i * GAP_PENALTY
        traceback[i][0] = 1

    for j in range(1, n + 1):
        dp[0][j] = j * GAP_PENALTY
        traceback[0][j] = 2

    for i in range(1, m + 1):
        ref_ch = reference[i - 1]
        prev_row = dp[i - 1]
        curr_row = dp[i]
        tb_row = traceback[i]
        for j in range(1, n + 1):
            diag = prev_row[j - 1] + _score_pair(ref_ch, query[j - 1])
            up = prev_row[j] + GAP_PENALTY
            left = curr_row[j - 1] + GAP_PENALTY

            if diag >= up and diag >= left:
                curr_row[j] = diag
                tb_row[j] = 0
            elif up >= left:
                curr_row[j] = up
                tb_row[j] = 1
            else:
                curr_row[j] = left
                tb_row[j] = 2

    ref_aligned_chars: List[str] = []
    query_aligned_chars: List[str] = []

    i, j = m, n
    while i > 0 or j > 0:
        move = traceback[i][j]
        if move == 0:
            ref_aligned_chars.append(reference[i - 1])
            query_aligned_chars.append(query[j - 1])
            i -= 1
            j -= 1
        elif move == 1:
            ref_aligned_chars.append(reference[i - 1])
            query_aligned_chars.append("-")
            i -= 1
        else:
            ref_aligned_chars.append("-")
            query_aligned_chars.append(query[j - 1])
            j -= 1

    ref_aligned = "".join(reversed(ref_aligned_chars))
    query_aligned = "".join(reversed(query_aligned_chars))
    return dp[m][n], ref_aligned, query_aligned


def calculate_identity(ref_aligned: str, query_aligned: str) -> float:
    total_cols = len(ref_aligned)
    if total_cols == 0:
        return 0.0
    exact_acgt_matches = sum(
        1
        for r, q in zip(ref_aligned, query_aligned)
        if r == q and r in {"A", "C", "G", "T"}
    )
    return round(exact_acgt_matches / total_cols, 6)


def extract_variants(ref_aligned: str, query_aligned: str) -> List[Dict]:
    variants: List[Dict] = []
    ref_pos = 0
    idx = 0
    total_cols = len(ref_aligned)

    while idx < total_cols:
        r_ch = ref_aligned[idx]
        q_ch = query_aligned[idx]

        if r_ch != "-" and q_ch != "-":
            ref_pos += 1
            if r_ch in {"A", "C", "G", "T"} and q_ch in {"A", "C", "G", "T"} and r_ch != q_ch:
                variants.append(
                    {
                        "position": ref_pos,
                        "ref": r_ch,
                        "alt": q_ch,
                        "type": "SNV",
                    }
                )
            idx += 1
        elif r_ch != "-" and q_ch == "-":
            start_pos = ref_pos + 1
            deleted_bases: List[str] = []
            while idx < total_cols and ref_aligned[idx] != "-" and query_aligned[idx] == "-":
                deleted_bases.append(ref_aligned[idx])
                ref_pos += 1
                idx += 1
            variants.append(
                {
                    "position": start_pos,
                    "ref": "".join(deleted_bases),
                    "alt": "-",
                    "type": "deletion",
                }
            )
        else:
            preceding_pos = ref_pos
            inserted_bases: List[str] = []
            while idx < total_cols and ref_aligned[idx] == "-" and query_aligned[idx] != "-":
                inserted_bases.append(query_aligned[idx])
                idx += 1
            variants.append(
                {
                    "position": preceding_pos,
                    "ref": "-",
                    "alt": "".join(inserted_bases),
                    "type": "insertion",
                }
            )

    variants.sort(key=lambda item: item["position"])
    return variants