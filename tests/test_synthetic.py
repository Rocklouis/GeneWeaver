import pytest

from pathlib import Path

from Bio import SeqIO

from geneweaver.synthetic import write_synthetic_fasta


def test_synthetic_fasta_length(tmp_path):
    output = tmp_path / "genome.fasta"

    written = write_synthetic_fasta(
        output,
        genome_length=10_000,
        seed=42,
    )

    record = next(SeqIO.parse(output, "fasta"))

    assert written == 10_000
    assert len(record.seq) == 10_000
    assert set(str(record.seq)) <= set("ACGT")


def test_synthetic_fasta_reproducibility(tmp_path):
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"

    write_synthetic_fasta(first, 5_000, seed=42)
    write_synthetic_fasta(second, 5_000, seed=42)

    first_record = next(SeqIO.parse(first, "fasta"))
    second_record = next(SeqIO.parse(second, "fasta"))

    assert str(first_record.seq) == str(second_record.seq)


def test_different_seeds_produce_different_sequences(tmp_path):
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"

    write_synthetic_fasta(first, 5_000, seed=42)
    write_synthetic_fasta(second, 5_000, seed=99)

    first_record = next(SeqIO.parse(first, "fasta"))
    second_record = next(SeqIO.parse(second, "fasta"))

    assert str(first_record.seq) != str(second_record.seq)


def test_invalid_genome_length(tmp_path):
    with pytest.raises(ValueError):
        write_synthetic_fasta(
            tmp_path / "invalid.fasta",
            genome_length=0,
        )


def test_invalid_block_size(tmp_path):
    with pytest.raises(ValueError):
        write_synthetic_fasta(
            tmp_path / "invalid.fasta",
            genome_length=100,
            block_size=0,
        )


def test_invalid_line_width(tmp_path):
    with pytest.raises(ValueError):
        write_synthetic_fasta(
            tmp_path / "invalid.fasta",
            genome_length=100,
            line_width=0,
        )