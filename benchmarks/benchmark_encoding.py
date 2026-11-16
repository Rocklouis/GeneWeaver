
import argparse
import time

from geneweaver.fasta import (
    encode_dna,
    iter_fasta_chunks,
)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
    )

    parser.add_argument(
        "--fasta",
        default="data/mock_genome.fasta",
    )

    args = parser.parse_args()

    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")

    total_bases = 0
    encoded_bases = 0

    start = time.perf_counter()

    for chunk in iter_fasta_chunks(
        args.fasta,
        chunk_size=args.chunk_size,
    ):
        encoded = encode_dna(chunk.sequence)

        total_bases += chunk.core_length
        encoded_bases += len(encoded)

        assert encoded.dtype.name == "uint8"

    elapsed = time.perf_counter() - start

    print("=" * 45)
    print("GENEWEAVER DNA ENCODING BENCHMARK")
    print("=" * 45)

    print(f"Original bases: {total_bases:,}")
    print(f"Encoded bases:  {encoded_bases:,}")
    print(f"Execution time: {elapsed:.6f} seconds")

    if elapsed > 0:
        print(
            f"Encoding throughput: "
            f"{encoded_bases / elapsed:,.0f} bases/second"
        )

    print("DNA ENCODING BENCHMARK COMPLETED")


if __name__ == "__main__":
    main()
