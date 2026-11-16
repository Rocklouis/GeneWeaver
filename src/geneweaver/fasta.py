
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO


DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4,
}

DNA_DECODING = "ACGTN"


@dataclass(frozen=True)
class GenomeChunk:
    record_id: str
    start: int
    sequence: str
    core_length: int

    @property
    def core_end(self) -> int:
        """Exclusive end coordinate of this chunk's core region."""
        return self.start + self.core_length


def normalize_dna(sequence: str) -> str:
    """Normalize DNA and replace unsupported bases with N."""
    sequence = sequence.upper().replace("U", "T")

    return "".join(
        base if base in DNA_ENCODING else "N"
        for base in sequence
    )


def validate_dna(sequence: str) -> bool:
    """Check whether all characters are valid DNA symbols."""
    return all(base in DNA_ENCODING for base in sequence.upper())


def encode_dna(sequence: str) -> np.ndarray:
    """
    Convert DNA into a NumPy uint8 array.

    A = 0, C = 1, G = 2, T = 3, N = 4
    """
    sequence = sequence.upper().replace("U", "T")

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    return np.fromiter(
        (DNA_ENCODING[base] for base in sequence),
        dtype=np.uint8,
        count=len(sequence),
    )


def decode_dna(encoded: np.ndarray) -> str:
    """Convert numerical DNA data back into a DNA string."""
    if np.any(encoded > 4):
        raise ValueError("Invalid encoded DNA value.")

    return "".join(
        DNA_DECODING[int(value)]
        for value in encoded
    )


def required_overlap(
    guide_length: int,
    pam_length: int = 3,
) -> int:
    """
    Calculate overlap needed to preserve complete alignment windows.

    A complete alignment window contains guide + PAM.
    """
    if guide_length <= 0 or pam_length <= 0:
        raise ValueError("Guide and PAM lengths must be positive.")

    return guide_length + pam_length - 1


def iter_fasta_chunks(
    fasta_path: str | Path,
    chunk_size: int = 1_000_000,
    overlap: int = 0,
):
    """
    Parse FASTA records and yield manageable genome chunks.

    Each chunk has a non-overlapping core region and may
    contain extra sequence at its end for alignment overlap.
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.is_file():
        raise FileNotFoundError(
            f"FASTA file not found: {fasta_path}"
        )

    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive.")

    if overlap < 0:
        raise ValueError("Overlap cannot be negative.")

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        sequence = normalize_dna(str(record.seq))

        for start in range(0, len(sequence), chunk_size):
            core_length = min(
                chunk_size,
                len(sequence) - start,
            )

            chunk_sequence = sequence[
                start:start + chunk_size + overlap
            ]

            yield GenomeChunk(
                record_id=record.id,
                start=start,
                sequence=chunk_sequence,
                core_length=core_length,
            )


def iter_alignment_chunks(
    fasta_path: str | Path,
    guide_length: int,
    pam_length: int = 3,
    chunk_size: int = 1_000_000,
):
    """
    Yield FASTA chunks with sufficient overlap for alignment.

    The core region determines which genomic start coordinates
    belong to each chunk. Downstream alignment code should use
    this region to avoid reporting duplicate hits.
    """
    overlap = required_overlap(
        guide_length,
        pam_length,
    )

    yield from iter_fasta_chunks(
        fasta_path,
        chunk_size=chunk_size,
        overlap=overlap,
    )
