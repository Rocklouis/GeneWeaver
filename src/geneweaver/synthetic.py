
from pathlib import Path
import random


DNA_BASES = "ACGT"


def write_synthetic_fasta(
    output_path: str | Path,
    genome_length: int,
    seed: int = 42,
    record_id: str = "synthetic_chromosome",
    block_size: int = 1_000_000,
    line_width: int = 80,
) -> int:
    """
    Generate a reproducible synthetic DNA FASTA file.

    DNA is generated in blocks to avoid holding the entire
    genome in memory.

    Args:
        output_path: Destination FASTA file.
        genome_length: Number of DNA bases to generate.
        seed: Random seed for reproducibility.
        record_id: FASTA record identifier.
        block_size: Number of bases generated per block.
        line_width: Number of bases per FASTA output line.

    Returns:
        Number of DNA bases written.
    """
    if genome_length <= 0:
        raise ValueError("Genome length must be positive.")

    if block_size <= 0:
        raise ValueError("Block size must be positive.")

    if line_width <= 0:
        raise ValueError("Line width must be positive.")

    if not record_id or any(c.isspace() for c in record_id):
        raise ValueError("Record ID must be non-empty and contain no spaces.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)

    with output_path.open("w", encoding="ascii", newline="\n") as fasta:
        fasta.write(f">{record_id}\n")

        remaining = genome_length

        while remaining > 0:
            current_block = min(block_size, remaining)

            sequence = "".join(
                DNA_BASES[rng.getrandbits(2)]
                for _ in range(current_block)
            )

            for start in range(0, len(sequence), line_width):
                fasta.write(
                    sequence[start:start + line_width] + "\n"
                )

            remaining -= current_block

    return genome_length
