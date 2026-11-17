
from Bio import SeqIO
import pytest

from geneweaver.fasta import iter_fasta_chunks
from geneweaver.streaming import iter_fasta_chunks_streaming


def test_streaming_matches_existing_parser(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(
        ">chr1\n"
        "ACGTACGTACGTACGTACGT\n"
        ">chr2\n"
        "TTTTCCCCAAAAGGGG\n"
    )

    expected = list(
        iter_fasta_chunks(fasta, chunk_size=8, overlap=3)
    )

    actual = list(
        iter_fasta_chunks_streaming(
            fasta,
            chunk_size=8,
            overlap=3,
        )
    )

    assert actual == expected


def test_streaming_preserves_coordinates(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1\nACGTACGTACGTACGT\n")

    chunks = list(
        iter_fasta_chunks_streaming(
            fasta,
            chunk_size=6,
            overlap=2,
        )
    )

    assert [chunk.start for chunk in chunks] == [0, 6, 12]
    assert [chunk.core_length for chunk in chunks] == [6, 6, 4]


def test_streaming_preserves_total_sequence(tmp_path):
    fasta = tmp_path / "genome.fasta"
    sequence = "ACGT" * 25
    fasta.write_text(f">chr1\n{sequence}\n")

    chunks = list(
        iter_fasta_chunks_streaming(
            fasta,
            chunk_size=17,
            overlap=5,
        )
    )

    reconstructed = "".join(
        chunk.sequence[:chunk.core_length]
        for chunk in chunks
    )

    assert reconstructed == sequence


def test_streaming_handles_empty_record(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text(">chr1\n")

    chunks = list(iter_fasta_chunks_streaming(fasta))

    assert chunks == []


def test_streaming_rejects_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        list(
            iter_fasta_chunks_streaming(
                tmp_path / "missing.fasta"
            )
        )


def test_streaming_rejects_invalid_chunk_size(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1\nACGT\n")

    with pytest.raises(ValueError):
        list(
            iter_fasta_chunks_streaming(
                fasta,
                chunk_size=0,
            )
        )


def test_streaming_rejects_invalid_overlap(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1\nACGT\n")

    with pytest.raises(ValueError):
        list(
            iter_fasta_chunks_streaming(
                fasta,
                overlap=-1,
            )
        )


def test_streaming_handles_sequence_without_header(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ACGTACGT\n")

    with pytest.raises(ValueError):
        list(iter_fasta_chunks_streaming(fasta))
