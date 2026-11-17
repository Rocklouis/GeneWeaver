
from pathlib import Path

from geneweaver.synthetic import write_synthetic_fasta


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "benchmarks"

DATASETS = [
    ("genome_100kb.fasta", 100_000, 42),
    ("genome_1mb.fasta", 1_000_000, 43),
    ("genome_10mb.fasta", 10_000_000, 44),
]


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 50)
    print("GENEWEAVER SYNTHETIC GENOME GENERATOR")
    print("=" * 50)

    for filename, genome_length, seed in DATASETS:
        output_path = DATA_DIR / filename

        written = write_synthetic_fasta(
            output_path=output_path,
            genome_length=genome_length,
            seed=seed,
            block_size=100_000,
        )

        file_size = output_path.stat().st_size

        print(f"\nDataset: {filename}")
        print(f"Bases generated: {written:,}")
        print(f"File size: {file_size:,} bytes")
        print("Status: SUCCESS")

    print("\nAll synthetic datasets generated.")


if __name__ == "__main__":
    main()
