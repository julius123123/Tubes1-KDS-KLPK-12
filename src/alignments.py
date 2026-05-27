from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class AlignmentResult:
    score: int
    aligned_a: str
    aligned_b: str

def _init_matrix(rows: int, cols: int, value: int = 0) -> List[List[int]]:
    return [[value for _ in range(cols)] for _ in range(rows)]

def needleman_wunsch(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2,) -> Tuple[int, str, str]:
    rows = len(a) + 1
    cols = len(b) + 1
    score = _init_matrix(rows, cols)
    trace = _init_matrix(rows, cols)

    for i in range(1, rows):
        score[i][0] = i * gap
        trace[i][0] = 1
    for j in range(1, cols):
        score[0][j] = j * gap
        trace[0][j] = 2

    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                trace[i][j] = 0
            elif best == up:
                trace[i][j] = 1
            else:
                trace[i][j] = 2

    i, j = len(a), len(b)
    aligned_a = []
    aligned_b = []
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 0:
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif direction == 1:
            aligned_a.append(a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(b[j - 1])
            j -= 1

    return score[-1][-1], "".join(reversed(aligned_a)), "".join(reversed(aligned_b))

def smith_waterman(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2, ) -> Tuple[int, str, str]:
    rows = len(a) + 1
    cols = len(b) + 1
    score = _init_matrix(rows, cols)
    trace = _init_matrix(rows, cols, -1)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(0, diag, up, left)
            score[i][j] = best
            if best == 0:
                trace[i][j] = -1
            elif best == diag:
                trace[i][j] = 0
            elif best == up:
                trace[i][j] = 1
            else:
                trace[i][j] = 2
            if best > max_score:
                max_score = best
                max_pos = (i, j)

    i, j = max_pos
    aligned_a = []
    aligned_b = []
    while i > 0 and j > 0 and score[i][j] != 0:
        direction = trace[i][j]
        if direction == 0:
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif direction == 1:
            aligned_a.append(a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(b[j - 1])
            j -= 1

    return max_score, "".join(reversed(aligned_a)), "".join(reversed(aligned_b))
