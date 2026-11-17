
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


def _count_mismatches(
    sequence: str,
    start: int,
    guide: str,
    max_mismatches: int,
) -> int:
    """Count mismatches and stop early when the limit is exceeded."""
    mismatches = 0

    for index, expected in enumerate(guide):
        if sequence[start + index] != expected:
            mismatches += 1

            if mismatches > max_mismatches:
                return mismatches

    return mismatches


def _scan_ngg(
    sequence: str,
    guide: str,
    max_mismatches: int,
    strand: str,
    reference_length: int,
) -> list[OffTargetHit]:
    """
    Efficiently scan standard NGG PAM sites.

    Search for GG first, then inspect only possible PAM positions.
    """
    hits = []
    guide_length = len(guide)
    window_length = guide_length + 3

    gg_position = sequence.find("GG")

    while gg_position != -1:
        pam_start = gg_position - 1
        start = pam_start - guide_length

        if (
            pam_start >= 0
            and start >= 0
            and gg_position + 2 <= len(sequence)
            and sequence[pam_start] in "ACGT"
            and start + window_length <= len(sequence)
        ):
            mismatches = _count_mismatches(
                sequence,
                start,
                guide,
                max_mismatches,
            )

            if mismatches <= max_mismatches:
                protospacer = sequence[start:start + guide_length]
                pam = sequence[pam_start:pam_start + 3]

                if strand == "+":
                    genomic_position = start
                else:
                    genomic_position = (
                        reference_length - start - window_length
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

        gg_position = sequence.find("GG", gg_position + 1)

    return hits


def _scan_generic_pam(
    sequence: str,
    guide: str,
    max_mismatches: int,
    pam_pattern: str,
    strand: str,
    reference_length: int,
) -> list[OffTargetHit]:
    """General scanning path for PAM patterns other than NGG."""
    hits = []
    guide_length = len(guide)
    pam_length = len(pam_pattern)
    window_length = guide_length + pam_length

    for start in range(len(sequence) - window_length + 1):
        pam_start = start + guide_length
        pam = sequence[pam_start:pam_start + pam_length]

        if not pam_matches(pam, pam_pattern):
            continue

        mismatches = _count_mismatches(
            sequence,
            start,
            guide,
            max_mismatches,
        )

        if mismatches <= max_mismatches:
            if strand == "+":
                genomic_position = start
            else:
                genomic_position = (
                    reference_length - start - window_length
                )

            hits.append(
                OffTargetHit(
                    position=genomic_position,
                    strand=strand,
                    mismatches=mismatches,
                    protospacer=sequence[start:start + guide_length],
                    pam=pam,
                )
            )

    return hits


def find_off_targets(
    reference: str,
    guide: str,
    max_mismatches: int = 3,
    pam_pattern: str = "NGG",
) -> list[OffTargetHit]:
    """
    Search both DNA strands for potential CRISPR off-target sites.

    Positions are zero-based coordinates on the original reference.
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

    reverse_reference = reverse_complement(reference)

    hits = []

    for strand, sequence in (
        ("+", reference),
        ("-", reverse_reference),
    ):
        if pam_pattern == "NGG":
            hits.extend(
                _scan_ngg(
                    sequence,
                    guide,
                    max_mismatches,
                    strand,
                    len(reference),
                )
            )
        else:
            hits.extend(
                _scan_generic_pam(
                    sequence,
                    guide,
                    max_mismatches,
                    pam_pattern,
                    strand,
                    len(reference),
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
