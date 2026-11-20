
import numpy as np
import pytest

from geneweaver.gpu.encoding import encode_dna_gpu


def test_encode_all_dna_bases():
    result = encode_dna_gpu("ACGTN")

    expected = np.array([0, 1, 2, 3, 4], dtype=np.uint8)

    np.testing.assert_array_equal(result, expected)


def test_encode_lowercase_dna():
    result = encode_dna_gpu("acgtn")

    expected = np.array([0, 1, 2, 3, 4], dtype=np.uint8)

    np.testing.assert_array_equal(result, expected)


def test_encode_mixed_case_dna():
    result = encode_dna_gpu("AcGtN")

    expected = np.array([0, 1, 2, 3, 4], dtype=np.uint8)

    np.testing.assert_array_equal(result, expected)


def test_invalid_base_encoding():
    result = encode_dna_gpu("ACXGT")

    expected = np.array([0, 1, 255, 2, 3], dtype=np.uint8)

    np.testing.assert_array_equal(result, expected)


def test_empty_sequence():
    result = encode_dna_gpu("")

    assert result.size == 0
    assert result.dtype == np.uint8


def test_large_dna_sequence():
    import random

    random.seed(42)

    sequence = "".join(
        random.choices("ACGTN", k=100_000)
    )

    mapping = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
        "N": 4,
    }

    expected = np.array(
        [mapping[base] for base in sequence],
        dtype=np.uint8,
    )

    result = encode_dna_gpu(sequence)

    np.testing.assert_array_equal(result, expected)


def test_invalid_input_type():
    with pytest.raises(TypeError, match="must be a string"):
        encode_dna_gpu(123)


def test_non_ascii_input():
    with pytest.raises(UnicodeEncodeError):
        encode_dna_gpu("ACGTé")
