from pathlib import Path

import matplotlib

from src.alignments import needleman_wunsch, smith_waterman
from src.utils import load_fasta
from src.classifier import classify_queries

DATA_DIR = Path("data")

def run_demo() -> None:
    matplotlib.use("Agg")
    refs = load_fasta(DATA_DIR / "refs.fasta")
    queries = load_fasta(DATA_DIR / "sequences.fasta")

    print("== Needleman-Wunsch (Global) Demo ==")
    s1_id, s1 = next(iter(queries.items()))
    s2_id, s2 = next(iter(refs.items()))
    score, a1, a2 = needleman_wunsch(s1, s2)
    print(f"Query: {s1_id}")
    print(f"Ref:   {s2_id}")
    print(f"Score: {score}")
    print(a1)
    print(a2)
    print()

    print("== Smith-Waterman (Local) Demo ==")
    score, a1, a2 = smith_waterman(s1, s2)
    print(f"Query: {s1_id}")
    print(f"Ref:   {s2_id}")
    print(f"Score: {score}")
    print(a1)
    print(a2)
    print()

    print("== Klasifikasi Varian ==")
    results = classify_queries(queries, refs)
    for q_id, result in results.items():
        print(f"{q_id} -> {result.label} (score={result.score}, sim={result.similarity:.3f})")
    print()


if __name__ == "__main__":
    run_demo()