
import pytest

from geneweaver.alignment import (
    find_off_targets,
    reverse_complement,
)


GUIDE = "ACGTACGTACGTACGTACGT"


def test_exact_match():
    reference = "TTT" + GUIDE + "AGG" + "TTT"

    hits = find_off_targets(reference, GUIDE, 0)

    assert any(
        hit.position == 3
        and hit.strand == "+"
        and hit.mismatches == 0
        for hit in hits
    )


def test_one_mismatch():
    mutated = GUIDE[:5] + "A" + GUIDE[6:]
    reference = "AAA" + mutated + "TGG" + "AAA"

    hits = find_off_targets(reference, GUIDE, 1)

    assert any(
        hit.position == 3 and hit.mismatches == 1
        for hit in hits
    )


def test_reverse_strand():
    target = GUIDE + "AGG"
    reverse = reverse_complement(target)

    reference = "TTT" + reverse + "TTT"

    hits = find_off_targets(reference, GUIDE, 0)

    assert any(
        hit.position == 3 and hit.strand == "-"
        for hit in hits
    )


def test_invalid_guide():
    with pytest.raises(ValueError):
        find_off_targets("ACGTACGT", "ACGTX")


def test_negative_mismatches():
    with pytest.raises(ValueError):
        find_off_targets("ACGTACGT", "ACGT", -1)


def test_wrong_pam():
    reference = "TTT" + GUIDE + "AAA" + "TTT"

    hits = find_off_targets(reference, GUIDE, 0)

    assert not hits
