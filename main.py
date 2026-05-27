from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo

from src.alignments import needleman_wunsch, smith_waterman
from src.classifier import classify_queries
from src.phylo import build_tree
from src.utils import load_fasta

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

    print("== Pohon Filogenetik (ASCII) ==")
    tree = build_tree({**refs, **queries})
    print(tree)
    out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)

    tree_path = out_dir / "tree.nwk"
    tree_path.write_text(tree.format("newick"))
    print(f"Newick saved to: {tree_path}")

    fig = plt.figure(figsize=(8, 6), dpi=160)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    png_path = out_dir / "tree.png"
    svg_path = out_dir / "tree.svg"
    fig.tight_layout()
    fig.savefig(png_path)
    fig.savefig(svg_path)
    print(f"PNG saved to: {png_path}")
    print(f"SVG saved to: {svg_path}")


if __name__ == "__main__":
    run_demo()
