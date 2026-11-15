
import pytest

from geneweaver.fasta import iter_fasta_chunks


def test_fasta_parsing_and_chunking(tmp_path):
    fasta = tmp_path / "mock.fasta"

    fasta.write_text(
        ">chr1\n"
        "ACGTAC\n"
        "GTNN\n"
        ">chr2\n"
        "TTTT\n"
    )

    chunks = list(
        iter_fasta_chunks(
            fasta,
            chunk_size=5,
            overlap=2,
        )
    )

    assert len(chunks) == 3

    assert chunks[0].record_id == "chr1"
    assert chunks[0].start == 0
    assert chunks[0].sequence == "ACGTACG"

    assert chunks[1].record_id == "chr1"
    assert chunks[1].start == 5
    assert chunks[1].sequence == "CGTNN"

    assert chunks[2].record_id == "chr2"
    assert chunks[2].sequence == "TTTT"


def test_ambiguous_bases_are_normalized(tmp_path):
    fasta = tmp_path / "ambiguous.fasta"
    fasta.write_text(">chr1\nACGTURYS\n")

    chunks = list(iter_fasta_chunks(fasta, chunk_size=20))

    assert chunks[0].sequence == "ACGTTNNN"


def test_missing_fasta():
    with pytest.raises(FileNotFoundError):
        list(iter_fasta_chunks("missing.fasta"))


def test_invalid_chunk_size(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">chr1\nACGT\n")

    with pytest.raises(ValueError):
        list(iter_fasta_chunks(fasta, chunk_size=0))
