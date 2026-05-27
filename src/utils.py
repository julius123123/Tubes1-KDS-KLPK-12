from pathlib import Path
from typing import Dict

from Bio import SeqIO


def load_fasta(path: Path) -> Dict[str, str]:
    records = SeqIO.parse(str(path), "fasta")
    result: Dict[str, str] = {}
    for record in records:
        result[record.id] = str(record.seq).upper()
    return result
