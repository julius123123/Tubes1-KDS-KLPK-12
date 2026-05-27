from typing import Dict, List

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from src.alignments import needleman_wunsch, smith_waterman


def _pairwise_distance(a: str, b: str, match: int, mismatch: int, gap: int) -> float:
    score, _, _ = needleman_wunsch(a, b, match, mismatch, gap)
    max_len = max(len(a), len(b))
    if max_len == 0:
        return 1.0
    sim = score / float(match * max_len)
    sim = max(0.0, min(1.0, sim))
    return 1.0 - sim

def _pairwise_distance_sw(a: str, b: str, match: int, mismatch: int, gap: int) -> float:
    score, align_a, align_b = smith_waterman(a, b, match, mismatch, gap)
    local_len = max(len(align_a.replace('-', '')), len(align_b.replace('-', '')))
    
    if local_len == 0:
        return 1.0
    
    sim = score / float(match * local_len)
    sim = max(0.0, min(1.0, sim))
    return 1.0 - sim


def build_tree(
    sequences: Dict[str, str],
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
    method: str = "nw"
) -> Phylo.BaseTree.Tree:
    names: List[str] = list(sequences.keys())
    matrix: List[List[float]] = []
    for i, name_i in enumerate(names):
        row: List[float] = []
        for j in range(i):
            if method == "nw":
                dist = _pairwise_distance(
                    sequences[name_i],
                    sequences[names[j]],
                    match,
                    mismatch,
                    gap,
                )
            elif method == "sw":
                dist = _pairwise_distance_sw(
                    sequences[name_i],
                    sequences[names[j]],
                    match,
                    mismatch,
                    gap,
                )
            row.append(dist)
        row.append(0.0)
        matrix.append(row)

    dm = DistanceMatrix(names=names, matrix=matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree
