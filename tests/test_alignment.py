
import pytest

from geneweaver.alignment import (
    find_off_targets,
    reverse_complement,
)


GUIDE = "ACGTACGTACGTACGTACGT"


def test_exact_match():
    reference = "AAAAA" + GUIDE + "AGG" + "TTTTT"

    hits = find_off_targets(
        reference=reference,
        guide=GUIDE,
        max_mismatches=0,
    )

    assert len(hits) >= 1

    matching_hits = [
        hit for hit in hits
        if hit.position == 5
        and hit.strand == "+"
        and hit.mismatches == 0
        and hit.pam == "AGG"
    ]

    assert len(matching_hits) == 1


def test_one_mismatch():
    mutated_guide = "TCGTACGTACGTACGTACGT"
    reference = "AAAAA" + mutated_guide + "AGG" + "TTTTT"

    hits = find_off_targets(
        reference=reference,
        guide=GUIDE,
        max_mismatches=1,
    )

    assert any(
        hit.position == 5
        and hit.strand == "+"
        and hit.mismatches == 1
        for hit in hits
    )


def test_reverse_strand():
    reverse_guide = reverse_complement(GUIDE)

    # A reverse-strand NGG site appears as CCN
    # followed by the reverse complement of the guide.
    reference = (
        "AAAAA"
        + "CCT"
        + reverse_guide
        + "TTTTT"
    )

    hits = find_off_targets(
        reference=reference,
        guide=GUIDE,
        max_mismatches=0,
    )

    assert any(
        hit.strand == "-"
        and hit.mismatches == 0
        for hit in hits
    )

def test_invalid_guide():
    with pytest.raises(ValueError):
        find_off_targets(
            reference="ACGTACGTACGT",
            guide="ACGTXYZ",
        )


def test_negative_mismatches():
    with pytest.raises(ValueError):
        find_off_targets(
            reference="ACGTACGTACGT",
            guide="ACGT",
            max_mismatches=-1,
        )


def test_wrong_pam():
    reference = "AAAAA" + GUIDE + "AAA" + "TTTTT"

    hits = find_off_targets(
        reference=reference,
        guide=GUIDE,
        max_mismatches=0,
    )

    assert not any(
        hit.position == 5 and hit.strand == "+"
        for hit in hits
    )
