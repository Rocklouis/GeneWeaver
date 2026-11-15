
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True)
class GenomeChunk:
    record_id: str
    start: int
    sequence: str


def iter_fasta_chunks(
    fasta_path: str | Path,
    chunk_size: int = 1_000_000,
    overlap: int = 0,
):
    """
    Read a FASTA file and yield manageable DNA chunks.

    Chunk start positions are zero-based.
    Overlap adds extra sequence to each chunk without
    changing the distance between consecutive chunk starts.
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive.")

    if overlap < 0:
        raise ValueError("Overlap cannot be negative.")

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        sequence = str(record.seq).upper().replace("U", "T")

        # Convert ambiguous and unsupported symbols to N.
        sequence = "".join(
            base if base in "ACGTN" else "N"
            for base in sequence
        )

        for start in range(0, len(sequence), chunk_size):
            chunk = sequence[
                start:start + chunk_size + overlap
            ]

            yield GenomeChunk(
                record_id=record.id,
                start=start,
                sequence=chunk,
            )
