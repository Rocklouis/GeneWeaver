
import gc
import time
import tracemalloc
from pathlib import Path

from geneweaver.fasta import iter_fasta_chunks
from geneweaver.streaming import iter_fasta_chunks_streaming


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FASTA_FILE = PROJECT_ROOT / "data" / "benchmarks" / "genome_10mb.fasta"

CHUNK_SIZE = 1_000_000
OVERLAP = 23


def benchmark_reader(name, reader):
    gc.collect()
    tracemalloc.start()

    start_time = time.perf_counter()
    total_bases = 0
    chunk_count = 0

    for chunk in reader(
        FASTA_FILE,
        chunk_size=CHUNK_SIZE,
        overlap=OVERLAP,
    ):
        total_bases += chunk.core_length
        chunk_count += 1

    elapsed = time.perf_counter() - start_time
    current_memory, peak_memory = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "name": name,
        "chunks": chunk_count,
        "bases": total_bases,
        "seconds": elapsed,
        "peak_mb": peak_memory / (1024 * 1024),
    }


def main():
    if not FASTA_FILE.is_file():
        raise FileNotFoundError(FASTA_FILE)

    print("=" * 60)
    print("GENEWEAVER: FASTA READER BENCHMARK")
    print("=" * 60)
    print(f"Dataset: {FASTA_FILE.name}")
    print(f"Chunk size: {CHUNK_SIZE:,}")
    print(f"Overlap: {OVERLAP}")
    print()

    results = [
        benchmark_reader("BioPython parser", iter_fasta_chunks),
        benchmark_reader("Streaming reader", iter_fasta_chunks_streaming),
    ]

    for result in results:
        print(f"\n{result['name']}")
        print("-" * 40)
        print(f"Chunks processed: {result['chunks']}")
        print(f"Bases processed:  {result['bases']:,}")
        print(f"Execution time:   {result['seconds']:.4f} seconds")
        print(f"Peak traced RAM:  {result['peak_mb']:.2f} MB")

    if results[0]["bases"] != results[1]["bases"]:
        raise RuntimeError("Readers processed different base counts!")

    if results[0]["chunks"] != results[1]["chunks"]:
        raise RuntimeError("Readers produced different chunk counts!")

    print("\nCorrectness check: PASS")
    print("Benchmark completed.")


if __name__ == "__main__":
    main()
