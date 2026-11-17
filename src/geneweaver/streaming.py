
from pathlib import Path

from geneweaver.fasta import GenomeChunk, normalize_dna


def iter_fasta_chunks_streaming(
    fasta_path: str | Path,
    chunk_size: int = 1_000_000,
    overlap: int = 0,
):
    """
    Read FASTA records line by line and yield overlapping chunks.

    Uses a list of sequence fragments to reduce repeated string
    concatenation. Does not load an entire FASTA record into memory.

    Chunk coordinates are zero-based and relative to each record.
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

    record_id = None
    parts = []
    buffer_length = 0
    total_bases = 0
    next_start = 0

    def finish_record():
        nonlocal parts, buffer_length, next_start

        buffer = "".join(parts)

        while next_start < total_bases:
            core_length = min(
                chunk_size,
                total_bases - next_start,
            )

            yield GenomeChunk(
                record_id=record_id,
                start=next_start,
                sequence=buffer[:core_length + overlap],
                core_length=core_length,
            )

            buffer = buffer[core_length:]
            next_start += core_length

        parts = []
        buffer_length = 0

    with fasta_path.open("r", encoding="ascii") as fasta:
        for line in fasta:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if record_id is not None:
                    yield from finish_record()

                header = line[1:].strip()

                if not header:
                    raise ValueError("FASTA record ID cannot be empty.")

                record_id = header.split()[0]
                parts = []
                buffer_length = 0
                total_bases = 0
                next_start = 0
                continue

            if record_id is None:
                raise ValueError(
                    "FASTA sequence encountered before record header."
                )

            normalized = normalize_dna(line)

            parts.append(normalized)
            buffer_length += len(normalized)
            total_bases += len(normalized)

            target_length = chunk_size + overlap

            if buffer_length >= target_length:
                buffer = "".join(parts)

                while len(buffer) >= target_length:
                    yield GenomeChunk(
                        record_id=record_id,
                        start=next_start,
                        sequence=buffer[:target_length],
                        core_length=chunk_size,
                    )

                    buffer = buffer[chunk_size:]
                    next_start += chunk_size

                parts = [buffer] if buffer else []
                buffer_length = len(buffer)

        if record_id is not None:
            yield from finish_record()
