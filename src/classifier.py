from dataclasses import dataclass
from typing import Dict, List, Tuple

from src.alignments import needleman_wunsch, smith_waterman


@dataclass
class Classification:
    label: str
    score: int
    similarity: float


def _similarity(score: int, max_len: int, match: int) -> float:
    if max_len == 0:
        return 0.0
    return max(0.0, min(1.0, score / float(match * max_len)))


def classify_queries(
    queries: Dict[str, str],
    refs: Dict[str, str],
    methods: List[str] = None,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> Dict[str, Dict[str, Classification]]:
    if methods is None:
        methods = ["Needleman Wunsch", "Smith Waterman"]

    all_results: Dict[str, Dict[str, Classification]] = {}
    for method in methods:
        method_results: Dict[str, Classification] = {}
        for q_id, q_seq in queries.items():
            best_label = "unknown"
            best_score = -10**9
            best_sim = 0.0
            for r_id, r_seq in refs.items():
                if method == "Needleman Wunsch":
                    score, _, _ = needleman_wunsch(q_seq, r_seq, match, mismatch, gap)
                elif method == "Smith Waterman":
                    score, _, _ = smith_waterman(q_seq, r_seq, match, mismatch, gap)
                else:
                    raise ValueError(f"Unknown method: {method}")
                sim = _similarity(score, max(len(q_seq), len(r_seq)), match)
                if score > best_score:
                    best_score = score
                    best_sim = sim
                    best_label = r_id.split("|")[0]
            method_results[q_id] = Classification(best_label, best_score, best_sim)
        all_results[method] = method_results

    return all_results
