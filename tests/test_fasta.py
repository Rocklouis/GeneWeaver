
import numpy as np
import pytest

from geneweaver.fasta import (
    decode_dna,
    encode_dna,
    iter_alignment_chunks,
    iter_fasta_chunks,
    normalize_dna,
    required_overlap,
    validate_dna,
)


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
    assert chunks[0].core_length == 5

    assert chunks[1].record_id == "chr1"
    assert chunks[1].start == 5
    assert chunks[1].sequence == "CGTNN"
    assert chunks[1].core_length == 5

    assert chunks[2].record_id == "chr2"
    assert chunks[2].sequence == "TTTT"


def test_ambiguous_bases_are_normalized(tmp_path):
    fasta = tmp_path / "ambiguous.fasta"
    fasta.write_text(">chr1\nACGTURYS\n")

    chunks = list(
        iter_fasta_chunks(fasta, chunk_size=20)
    )

    assert chunks[0].sequence == "ACGTTNNN"


def test_missing_fasta():
    with pytest.raises(FileNotFoundError):
        list(iter_fasta_chunks("missing.fasta"))


def test_invalid_chunk_size(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">chr1\nACGT\n")

    with pytest.raises(ValueError):
        list(iter_fasta_chunks(fasta, chunk_size=0))


def test_dna_encoding():
    encoded = encode_dna("ACGTN")

    assert encoded.dtype == np.uint8
    assert encoded.tolist() == [0, 1, 2, 3, 4]


def test_dna_decoding():
    encoded = np.array(
        [0, 1, 2, 3, 4],
        dtype=np.uint8,
    )

    assert decode_dna(encoded) == "ACGTN"


def test_encoding_round_trip():
    original = "ACGTNACGT"
    encoded = encode_dna(original)

    assert decode_dna(encoded) == original


def test_invalid_encoding():
    with pytest.raises(ValueError):
        encode_dna("ACGTXYZ")


def test_dna_normalization():
    assert normalize_dna("acgturys") == "ACGTTNNN"


def test_dna_validation():
    assert validate_dna("ACGTN")
    assert not validate_dna("ACGTXYZ")


def test_required_overlap():
    assert required_overlap(20, 3) == 22


def test_invalid_overlap_parameters():
    with pytest.raises(ValueError):
        required_overlap(0, 3)


def test_alignment_chunks_overlap(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1\nACGTACGTACGTACGT\n")

    chunks = list(
        iter_alignment_chunks(
            fasta,
            guide_length=4,
            pam_length=3,
            chunk_size=5,
        )
    )

    assert len(chunks) == 4
    assert chunks[0].start == 0
    assert chunks[0].sequence == "ACGTACGTACG"
    assert chunks[0].core_length == 5
    assert chunks[0].core_end == 5

    assert chunks[1].start == 5
    assert chunks[1].core_length == 5
    assert chunks[1].core_end == 10


def test_empty_fasta(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    assert list(iter_fasta_chunks(fasta)) == []
