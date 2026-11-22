
"""Parallel CRISPR off-target alignment using an OpenCL GPU."""

from dataclasses import dataclass

import numpy as np
import pyopencl as cl


@dataclass(frozen=True)
class GPUOffTargetHit:
    position: int
    strand: str
    mismatches: int
    protospacer: str
    pam: str


ALIGNMENT_KERNEL = """
__kernel void scan_alignment(
    __global const uchar *reference,
    __global const uchar *guide,
    __global const uchar *pam_pattern,
    __global int *hit_flags,
    __global int *mismatch_counts,
    const int reference_length,
    const int guide_length,
    const int pam_length,
    const int max_mismatches
)
{
    int start = get_global_id(0);
    int window_length = guide_length + pam_length;

    if (start + window_length > reference_length)
        return;

    int mismatches = 0;

    // Compare the protospacer against the guide.
    for (int i = 0; i < guide_length; i++)
    {
        if (reference[start + i] != guide[i])
        {
            mismatches++;

            if (mismatches > max_mismatches)
            {
                hit_flags[start] = 0;
                mismatch_counts[start] = mismatches;
                return;
            }
        }
    }

    // Validate the PAM.
    int pam_start = start + guide_length;

    for (int i = 0; i < pam_length; i++)
    {
        uchar observed = reference[pam_start + i];
        uchar expected = pam_pattern[i];

        // Code 4 represents N: it matches only A, C, G, or T.
        if (expected == 4)
        {
            if (observed > 3)
            {
                hit_flags[start] = 0;
                mismatch_counts[start] = mismatches;
                return;
            }
        }
        else if (observed != expected)
        {
            hit_flags[start] = 0;
            mismatch_counts[start] = mismatches;
            return;
        }
    }

    hit_flags[start] = 1;
    mismatch_counts[start] = mismatches;
}
"""


_BASE_TO_CODE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4,
}


def _encode(sequence):
    """Encode DNA bases as uint8 values."""
    return np.array(
        [_BASE_TO_CODE[base] for base in sequence],
        dtype=np.uint8,
    )


def _get_opencl_gpu():
    """Find the first available OpenCL GPU."""
    for platform in cl.get_platforms():
        devices = platform.get_devices(
            device_type=cl.device_type.GPU
        )

        if devices:
            return devices[0]

    raise RuntimeError("No OpenCL GPU device available.")


def _scan_one_strand(
    sequence,
    guide,
    pam_pattern,
    max_mismatches,
):
    """Scan all possible alignment positions on one strand."""

    guide_length = len(guide)
    pam_length = len(pam_pattern)
    window_length = guide_length + pam_length
    reference_length = len(sequence)

    if reference_length < window_length:
        return []

    reference_data = _encode(sequence)
    guide_data = _encode(guide)
    pam_data = _encode(pam_pattern)

    device = _get_opencl_gpu()
    context = cl.Context(devices=[device])
    queue = cl.CommandQueue(context)

    program = cl.Program(
        context,
        ALIGNMENT_KERNEL,
    ).build()

    candidate_count = reference_length - window_length + 1

    hit_flags = np.zeros(candidate_count, dtype=np.int32)
    mismatch_counts = np.zeros(candidate_count, dtype=np.int32)

    flags_buffer = cl.Buffer(
        context,
        cl.mem_flags.WRITE_ONLY,
        hit_flags.nbytes,
    )

    mismatch_buffer = cl.Buffer(
        context,
        cl.mem_flags.WRITE_ONLY,
        mismatch_counts.nbytes,
    )

    reference_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=reference_data,
    )

    guide_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=guide_data,
    )

    pam_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=pam_data,
    )

    program.scan_alignment(
        queue,
        (candidate_count,),
        None,
        reference_buffer,
        guide_buffer,
        pam_buffer,
        flags_buffer,
        mismatch_buffer,
        np.int32(reference_length),
        np.int32(guide_length),
        np.int32(pam_length),
        np.int32(max_mismatches),
    )

    cl.enqueue_copy(queue, hit_flags, flags_buffer)
    cl.enqueue_copy(queue, mismatch_counts, mismatch_buffer).wait()

    hits = []

    for start in np.flatnonzero(hit_flags):
        start = int(start)
        pam_start = start + guide_length

        hits.append(
            (
                start,
                int(mismatch_counts[start]),
                sequence[start:pam_start],
                sequence[pam_start:pam_start + pam_length],
            )
        )

    return hits


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    translation = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(translation)[::-1]


def find_off_targets_gpu(
    reference,
    guide,
    max_mismatches=3,
    pam_pattern="NGG",
):
    """
    Find CRISPR off-target sites using parallel GPU alignment.

    Returns GPUOffTargetHit objects sorted by original reference
    position, strand, and mismatch count.
    """
    if not isinstance(reference, str) or not isinstance(guide, str):
        raise TypeError("Reference and guide must be strings.")

    if not isinstance(pam_pattern, str):
        raise TypeError("PAM pattern must be a string.")

    reference = reference.upper().replace("U", "T")
    guide = guide.upper().replace("U", "T")
    pam_pattern = pam_pattern.upper()

    if not guide or any(base not in "ACGT" for base in guide):
        raise ValueError("Guide must contain only A, C, G and T.")

    if not pam_pattern or any(
        base not in "ACGTN" for base in pam_pattern
    ):
        raise ValueError("Invalid PAM pattern.")

    if not isinstance(max_mismatches, int) or isinstance(
        max_mismatches, bool
    ):
        raise TypeError("Maximum mismatches must be an integer.")

    if max_mismatches < 0:
        raise ValueError("Maximum mismatches cannot be negative.")

    if any(base not in "ACGTN" for base in reference):
        raise ValueError("Reference contains invalid DNA characters.")

    reverse_reference = reverse_complement(reference)
    reference_length = len(reference)
    window_length = len(guide) + len(pam_pattern)

    hits = []

    for strand, sequence in (
        ("+", reference),
        ("-", reverse_reference),
    ):
        strand_hits = _scan_one_strand(
            sequence,
            guide,
            pam_pattern,
            max_mismatches,
        )

        for start, mismatches, protospacer, pam in strand_hits:
            if strand == "+":
                genomic_position = start
            else:
                genomic_position = (
                    reference_length - start - window_length
                )

            hits.append(
                GPUOffTargetHit(
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
