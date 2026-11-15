
from dataclasses import dataclass


@dataclass(frozen=True)
class OffTargetHit:
    position: int
    strand: str
    mismatches: int
    protospacer: str
    pam: str


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    translation = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(translation)[::-1]


def normalize_dna(sequence: str) -> str:
    """Normalize DNA input and convert U to T."""
    return sequence.upper().replace("U", "T")


def pam_matches(sequence: str, pattern: str = "NGG") -> bool:
    """Check whether a DNA sequence matches a PAM pattern."""
    if len(sequence) != len(pattern):
        return False

    for base, expected in zip(sequence, pattern):
        if expected == "N":
            if base not in "ACGT":
                return False
        elif base != expected:
            return False

    return True


def find_off_targets(
    reference: str,
    guide: str,
    max_mismatches: int = 3,
    pam_pattern: str = "NGG",
) -> list[OffTargetHit]:
    """
    Search both DNA strands for potential CRISPR off-target sites.

    Positions are zero-based coordinates on the original reference.
    The PAM is located immediately after the protospacer in the
    orientation being scanned.
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
