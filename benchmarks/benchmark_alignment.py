
import gc
import statistics
import time
from pathlib import Path

from Bio import SeqIO

from geneweaver.alignment import find_off_targets
from geneweaver.gpu.alignment import find_off_targets_gpu


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "benchmarks"

DATASETS = [
    ("100 KB", DATA_DIR / "genome_100kb.fasta"),
    ("1 MB", DATA_DIR / "genome_1mb.fasta"),
    ("10 MB", DATA_DIR / "genome_10mb.fasta"),
]

def choose_guide(reference):
    """Find a 20-base guide followed by a valid NGG PAM."""
    for position in range(len(reference) - 22):
        guide = reference[position:position + 20]
        pam = reference[position + 20:position + 23]

        if (
            all(base in "ACGT" for base in guide)
            and pam[0] in "ACGT"
            and pam[1:] == "GG"
        ):
            return guide

    raise ValueError("No valid NGG target found in reference.")

MAX_MISMATCHES = 2
PAM_PATTERN = "NGG"

CPU_REPEATS = 3
GPU_REPEATS = 3


def load_reference(path):
    """Load and normalize a FASTA reference."""
    records = list(SeqIO.parse(str(path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            f"Expected exactly one FASTA record in {path.name}"
        )

    return str(records[0].seq).upper().replace("U", "T")


def normalize_hits(hits):
    """Convert CPU/GPU hit objects into comparable tuples."""
    return [
        (
            hit.position,
            hit.strand,
            hit.mismatches,
            hit.protospacer,
            hit.pam,
        )
        for hit in hits
    ]


def measure(function, repeats):
    """Measure repeated execution and return elapsed times."""
    timings = []

    for _ in range(repeats):
        gc.collect()

        start = time.perf_counter()
        result = function()
        elapsed = time.perf_counter() - start

        timings.append(elapsed)

    return result, timings


def benchmark_dataset(label, path):
    print("\n" + "=" * 65)
    print(f"Dataset: {label}")
    print("=" * 65)

    reference = load_reference(path)
    guide = choose_guide(reference)

    print(f"Reference length: {len(reference):,} bases")
    print(f"Guide: {guide}")
    print(f"Maximum mismatches: {MAX_MISMATCHES}")
    print(f"PAM: {PAM_PATTERN}")

    cpu_function = lambda: find_off_targets(
        reference,
        guide,
        MAX_MISMATCHES,
        PAM_PATTERN,
    )

    gpu_function = lambda: find_off_targets_gpu(
        reference,
        guide,
        MAX_MISMATCHES,
        PAM_PATTERN,
    )

    # CPU benchmark
    cpu_hits, cpu_times = measure(
        cpu_function,
        CPU_REPEATS,
    )

    # GPU first call includes OpenCL setup and kernel compilation.
    gpu_first_start = time.perf_counter()
    gpu_first_hits = gpu_function()
    gpu_first_time = time.perf_counter() - gpu_first_start

    # Verify correctness before reporting performance.
    if normalize_hits(cpu_hits) != normalize_hits(gpu_first_hits):
        raise AssertionError(
            f"CPU/GPU results differ for {label}"
        )

    # Warm-up and repeated GPU execution.
    gpu_function()

    gpu_hits, gpu_times = measure(
        gpu_function,
        GPU_REPEATS,
    )

    if normalize_hits(cpu_hits) != normalize_hits(gpu_hits):
        raise AssertionError(
            f"CPU/GPU results differ after warm-up for {label}"
        )

    cpu_average = statistics.mean(cpu_times)
    gpu_average = statistics.mean(gpu_times)

    speedup = cpu_average / gpu_average

    print(f"\nCPU hits: {len(cpu_hits)}")
    print(f"GPU hits: {len(gpu_hits)}")
    print("Correctness: PASS")

    print("\nExecution times:")
    print(
        f"CPU average: {cpu_average:.6f} seconds"
    )
    print(
        f"GPU first call: {gpu_first_time:.6f} seconds"
    )
    print(
        f"GPU warm average: {gpu_average:.6f} seconds"
    )

    print(f"\nWarm GPU speedup: {speedup:.3f}x")

    if speedup > 1:
        print("Result: GPU is faster in this measurement.")
    elif speedup < 1:
        print("Result: CPU is faster in this measurement.")
    else:
        print("Result: Equal execution time.")

    return {
        "dataset": label,
        "bases": len(reference),
        "hits": len(cpu_hits),
        "cpu_average_seconds": cpu_average,
        "gpu_first_seconds": gpu_first_time,
        "gpu_warm_average_seconds": gpu_average,
        "warm_speedup": speedup,
    }


def main():
    print("GeneWeaver CPU vs GPU Alignment Benchmark")
    print("==========================================")

    results = []

    for label, path in DATASETS:
        if not path.is_file():
            raise FileNotFoundError(path)

        results.append(
            benchmark_dataset(label, path)
        )

    print("\n\nFINAL BENCHMARK SUMMARY")
    print("=" * 90)

    print(
        f"{'Dataset':<12}"
        f"{'Bases':>13}"
        f"{'CPU (s)':>14}"
        f"{'GPU warm (s)':>17}"
        f"{'Speedup':>13}"
    )

    for result in results:
        print(
            f"{result['dataset']:<12}"
            f"{result['bases']:>13,}"
            f"{result['cpu_average_seconds']:>14.6f}"
            f"{result['gpu_warm_average_seconds']:>17.6f}"
            f"{result['warm_speedup']:>12.3f}x"
        )

    print("\nAll datasets completed successfully.")


if __name__ == "__main__":
    main()
