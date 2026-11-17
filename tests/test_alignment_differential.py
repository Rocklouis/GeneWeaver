
import random

import pytest

from geneweaver.alignment import (
    OffTargetHit,
    find_off_targets,
    normalize_dna,
    pam_matches,
    reverse_complement,
)


def reference_find_off_targets(
    reference,
    guide,
    max_mismatches=3,
    pam_pattern="NGG",
):
    """
    Original straightforward alignment algorithm.

    This acts as the reference against which the optimized
    implementation is tested.
    """
    reference = normalize_dna(reference)
    guide = normalize_dna(guide)
    pam_pattern = pam_pattern.upper()

    if not guide or any(base not in "ACGT" for base in guide):
        raise ValueError("Guide must contain only A, C, G and T.")

    if not pam_pattern or any(
        base not in "ACGTN" for base in pam_pattern
    ):
        raise ValueError("Invalid PAM pattern.")

    if max_mismatches < 0:
        raise ValueError("Maximum mismatches cannot be negative.")

    if any(base not in "ACGTN" for base in reference):
        raise ValueError("Reference contains invalid DNA characters.")

    hits = []
    window_length = len(guide) + len(pam_pattern)

    strands = [
        ("+", reference),
        ("-", reverse_complement(reference)),
    ]

    for strand, sequence in strands:
        for start in range(len(sequence) - window_length + 1):
            protospacer = sequence[start:start + len(guide)]

            pam_start = start + len(guide)
            pam = sequence[
                pam_start:pam_start + len(pam_pattern)
            ]

            if not pam_matches(pam, pam_pattern):
                continue

            mismatches = sum(
                a != b for a, b in zip(protospacer, guide)
            )

            if mismatches <= max_mismatches:
                if strand == "+":
                    genomic_position = start
                else:
                    genomic_position = (
                        len(reference) - start - window_length
                    )

                hits.append(
                    OffTargetHit(
                        position=genomic_position,
                        strand=strand,
                        mismatches=mismatches,
                        protospacer=protospacer,
                        pam=pam,
                    )
                )

    return sorted(
        hits,
        key=lambda hit: (
            hit.position,
            hit.strand,
            hit.mismatches,
        ),
    )


def test_randomized_ngg_differential():
    """Compare both algorithms on 100 reproducible random genomes."""
    rng = random.Random(12345)
    guide = "ACGTACGTACGTACGTACGT"

    for case in range(100):
        reference = "".join(
            rng.choices("ACGT", k=rng.randint(50, 500))
        )

        mismatch_limit = rng.randint(0, 5)

        expected = reference_find_off_targets(
            reference,
            guide,
            max_mismatches=mismatch_limit,
            pam_pattern="NGG",
        )

        actual = find_off_targets(
            reference,
            guide,
            max_mismatches=mismatch_limit,
            pam_pattern="NGG",
        )

        assert actual == expected, (
            f"NGG differential mismatch in case {case}"
        )


@pytest.mark.parametrize(
    "pam_pattern",
    ["NGG", "NAG", "TTA"],
)
def test_custom_pam_differential(pam_pattern):
    """Verify standard and custom PAM patterns."""
    rng = random.Random(54321)
    guide = "ACGTACGTACGT"

    for case in range(30):
        reference = "".join(
            rng.choices("ACGT", k=200)
        )

        expected = reference_find_off_targets(
            reference,
            guide,
            max_mismatches=3,
            pam_pattern=pam_pattern,
        )

        actual = find_off_targets(
            reference,
            guide,
            max_mismatches=3,
            pam_pattern=pam_pattern,
        )

        assert actual == expected, (
            f"PAM {pam_pattern} differential mismatch "
            f"in case {case}"
        )


def test_randomized_sequences_with_ambiguous_bases():
    """Verify that ambiguous N bases behave identically."""
    rng = random.Random(98765)
    guide = "ACGTACGTACGT"

    for case in range(30):
        reference = "".join(
            rng.choices("ACGTN", k=250)
        )

        expected = reference_find_off_targets(
            reference,
            guide,
            max_mismatches=4,
        )

        actual = find_off_targets(
            reference,
            guide,
            max_mismatches=4,
        )

        assert actual == expected, (
            f"Ambiguous-base mismatch in case {case}"
        )


def test_inserted_exact_match_on_both_orientations():
    """Verify known exact matches on forward and reverse strands."""
    guide = "ACGTACGTACGT"
    reverse_guide = reverse_complement(guide)

    reference = (
        "AAAAA"
        + guide + "AGG"
        + "CCCCCC"
        + reverse_guide + "AGG"
        + "TTTTT"
    )

    expected = reference_find_off_targets(
        reference,
        guide,
        max_mismatches=0,
    )

    actual = find_off_targets(
        reference,
        guide,
        max_mismatches=0,
    )

    assert actual == expected
    assert {hit.strand for hit in actual} == {"+", "-"}


def test_rna_guide_normalization():
    """Verify that U in a guide is normalized to T."""
    reference = "AAAAAACGTACGTACGTAGGTTTTT"

    expected = reference_find_off_targets(
        reference,
        "ACGUACGUACGU",
        max_mismatches=0,
    )

    actual = find_off_targets(
        reference,
        "ACGUACGUACGU",
        max_mismatches=0,
    )

    assert actual == expected
