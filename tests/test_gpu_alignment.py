
import random

import pytest

from geneweaver.alignment import find_off_targets
from geneweaver.gpu.alignment import find_off_targets_gpu


def normalize_hits(hits):
    """Convert hits into comparable tuples."""
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


def assert_cpu_gpu_match(reference, guide, mismatches=0, pam="NGG"):
    cpu_hits = find_off_targets(
        reference,
        guide,
        mismatches,
        pam,
    )

    gpu_hits = find_off_targets_gpu(
        reference,
        guide,
        mismatches,
        pam,
    )

    assert normalize_hits(gpu_hits) == normalize_hits(cpu_hits)


def test_gpu_forward_exact_match():
    assert_cpu_gpu_match(
        "AAACCCGGGTTTACGTAGG",
        "ACGT",
    )

    hits = find_off_targets_gpu(
        "AAACCCGGGTTTACGTAGG",
        "ACGT",
        0,
    )

    assert any(hit.strand == "+" for hit in hits)


def test_gpu_reverse_strand_match():
    reference = "CCTACGT"

    assert_cpu_gpu_match(reference, "ACGT")

    hits = find_off_targets_gpu(reference, "ACGT", 0)

    assert len(hits) == 1
    assert hits[0].position == 0
    assert hits[0].strand == "-"
    assert hits[0].mismatches == 0


@pytest.mark.parametrize(
    "pam,reference",
    [
        ("NGG", "ACGTAGG"),
        ("NAG", "ACGTAAG"),
        ("TTA", "ACGTTTA"),
    ],
)
def test_gpu_custom_pam_patterns(pam, reference):
    assert_cpu_gpu_match(
        reference,
        "ACGT",
        0,
        pam,
    )

    hits = find_off_targets_gpu(
        reference,
        "ACGT",
        0,
        pam,
    )

    assert len(hits) == 1


def test_gpu_one_mismatch():
    reference = "ACGTTAGG"

    assert_cpu_gpu_match(
        reference,
        "ACGT",
        1,
        "NGG",
    )


def test_gpu_zero_mismatch_rejects_mismatch():
    reference = "ACGTTAGG"

    assert_cpu_gpu_match(
        reference,
        "ACGT",
        0,
        "NGG",
    )

    assert find_off_targets_gpu(
        reference,
        "ACGT",
        0,
        "NGG",
    ) == []


def test_gpu_empty_reference():
    assert_cpu_gpu_match("", "ACGT")


def test_gpu_short_reference():
    assert_cpu_gpu_match("ACG", "ACGT")


def test_gpu_invalid_guide():
    with pytest.raises(ValueError):
        find_off_targets_gpu("ACGTAGG", "ACGN")


def test_gpu_invalid_pam():
    with pytest.raises(ValueError):
        find_off_targets_gpu("ACGTAGG", "ACGT", 0, "XYZ")


def test_gpu_negative_mismatches():
    with pytest.raises(ValueError):
        find_off_targets_gpu("ACGTAGG", "ACGT", -1)


def test_gpu_invalid_reference():
    with pytest.raises(ValueError):
        find_off_targets_gpu("ACGXAGG", "ACGT")


def test_gpu_randomized_differential():
    random.seed(42)

    for _ in range(20):
        reference = "".join(
            random.choices("ACGT", k=100)
        )

        guide = "".join(
            random.choices("ACGT", k=12)
        )

        mismatches = random.randint(0, 3)

        assert_cpu_gpu_match(
            reference,
            guide,
            mismatches,
            "NGG",
        )
