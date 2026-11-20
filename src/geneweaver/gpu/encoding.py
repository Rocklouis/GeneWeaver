
"""GPU-accelerated DNA sequence encoding using OpenCL."""

import numpy as np
import pyopencl as cl


DNA_ENCODING_KERNEL = """
__kernel void encode_dna(
    __global const char *sequence,
    __global unsigned char *encoded,
    const int n
)
{
    int i = get_global_id(0);

    if (i < n)
    {
        char base = sequence[i];

        if (base == 'A' || base == 'a')
            encoded[i] = 0;
        else if (base == 'C' || base == 'c')
            encoded[i] = 1;
        else if (base == 'G' || base == 'g')
            encoded[i] = 2;
        else if (base == 'T' || base == 't')
            encoded[i] = 3;
        else if (base == 'N' || base == 'n')
            encoded[i] = 4;
        else
            encoded[i] = 255;
    }
}
"""


def get_opencl_gpu():
    """Return the first available OpenCL GPU."""
    for platform in cl.get_platforms():
        devices = platform.get_devices(
            device_type=cl.device_type.GPU
        )

        if devices:
            return devices[0]

    raise RuntimeError("No OpenCL GPU device available.")


def encode_dna_gpu(sequence):
    """
    Encode DNA bases on the GPU.

    Encoding:
        A = 0
        C = 1
        G = 2
        T = 3
        N = 4
        Invalid characters = 255

    The input is transferred to GPU memory, encoded in parallel,
    and the result is copied back to CPU memory.
    """
    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")

    if not sequence:
        return np.array([], dtype=np.uint8)

    sequence_bytes = np.frombuffer(
        sequence.encode("ascii"),
        dtype=np.uint8,
    )

    device = get_opencl_gpu()

    context = cl.Context(devices=[device])
    queue = cl.CommandQueue(context)

    program = cl.Program(
        context,
        DNA_ENCODING_KERNEL,
    ).build()

    encoded = np.empty(
        len(sequence_bytes),
        dtype=np.uint8,
    )

    sequence_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=sequence_bytes,
    )

    encoded_buffer = cl.Buffer(
        context,
        cl.mem_flags.WRITE_ONLY,
        encoded.nbytes,
    )

    program.encode_dna(
        queue,
        (len(sequence_bytes),),
        None,
        sequence_buffer,
        encoded_buffer,
        np.int32(len(sequence_bytes)),
    )

    cl.enqueue_copy(
        queue,
        encoded,
        encoded_buffer,
    ).wait()

    return encoded
