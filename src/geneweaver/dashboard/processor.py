
from pathlib import Path
from time import perf_counter
from typing import Callable

from geneweaver.alignment import find_off_targets
from geneweaver.fasta import required_overlap
from geneweaver.gpu.alignment import find_off_targets_gpu
from geneweaver.streaming import iter_fasta_chunks_streaming


def process_fasta_alignment(
    fasta_path: str | Path,
    guide: str,
    engine: str = "cpu",
    max_mismatches: int = 3,
    pam_pattern: str = "NGG",
    chunk_size: int = 250_000,
    progress_callback: Callable | None = None,
) -> dict:
    """Align a guide against a FASTA file using overlapping chunks."""

    guide = guide.strip().upper().replace("U", "T")
    pam_pattern = pam_pattern.strip().upper()

    overlap = required_overlap(
        len(guide),
        len(pam_pattern),
    )

    if engine == "cpu":
        aligner = find_off_targets
    elif engine == "gpu":
        aligner = find_off_targets_gpu
    else:
        raise ValueError("Engine must be 'cpu' or 'gpu'.")

    # First pass: calculate total work for progress reporting.
    total_chunks = 0
    total_bases = 0

    for chunk in iter_fasta_chunks_streaming(
        fasta_path,
        chunk_size=chunk_size,
        overlap=overlap,
    ):
        total_chunks += 1
        total_bases += chunk.core_length

    processed_chunks = 0
    processed_bases = 0
    all_hits = []

    start_time = perf_counter()

    # Second pass: perform actual alignment.
    for chunk in iter_fasta_chunks_streaming(
        fasta_path,
        chunk_size=chunk_size,
        overlap=overlap,
    ):
        chunk_hits = aligner(
            chunk.sequence,
            guide,
            max_mismatches,
            pam_pattern,
        )

        # Only count hits whose starting position belongs
        # to this chunk's non-overlapping core region.
        for hit in chunk_hits:
            if 0 <= hit.position < chunk.core_length:
                all_hits.append({
                    "record_id": chunk.record_id,
                    "position": chunk.start + hit.position,
                    "strand": hit.strand,
                    "mismatches": hit.mismatches,
                    "protospacer": hit.protospacer,
                    "pam": hit.pam,
                })

        processed_chunks += 1
        processed_bases += chunk.core_length

        if progress_callback is not None:
            progress_callback(
                processed_chunks,
                total_chunks,
                processed_bases,
                total_bases,
                len(all_hits),
            )

    elapsed = perf_counter() - start_time

    return {
        "engine": engine,
        "chunks_processed": processed_chunks,
        "bases_processed": processed_bases,
        "total_hits": len(all_hits),
        "hits": all_hits,
        "elapsed_seconds": elapsed,
    }
