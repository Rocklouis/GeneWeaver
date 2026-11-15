
import argparse
import random
import time

from geneweaver.alignment import find_off_targets


GUIDE = "ACGTACGTACGTACGTACGT"


def generate_genome(size: int) -> str:
    rng = random.Random(42)

    if size < len(GUIDE) + 3:
        raise ValueError("Genome size is too small.")

    genome = "".join(
        rng.choices("ACGT", k=size)
    )

    position = size // 2

    genome = (
        genome[:position]
        + GUIDE
        + "AGG"
        + genome[position + len(GUIDE) + 3:]
    )

    return genome


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bases",
        type=int,
        default=1_000_000,
    )

    args = parser.parse_args()

    if args.bases <= 0:
        parser.error("--bases must be positive")

    genome = generate_genome(args.bases)

    print("=" * 45)
    print("GENEWEAVER CPU BASELINE BENCHMARK")
    print("=" * 45)

    print(f"Genome size: {len(genome):,} base pairs")
    print(f"Guide length: {len(GUIDE)}")
    print("Maximum mismatches: 3")

    start = time.perf_counter()

    hits = find_off_targets(
        reference=genome,
        guide=GUIDE,
        max_mismatches=3,
    )

    end = time.perf_counter()

    elapsed = end - start

    print("\nRESULTS")
    print(f"Execution time: {elapsed:.6f} seconds")
    print(f"Total hits: {len(hits)}")

    if hits:
        print("\nFirst five candidate sites:")

        for hit in hits[:5]:
            print(
                f"Position: {hit.position}, "
                f"Strand: {hit.strand}, "
                f"Mismatches: {hit.mismatches}, "
                f"PAM: {hit.pam}"
            )

    print("\nCPU BASELINE COMPLETED")


if __name__ == "__main__":
    main()
