
"""Basic GPU kernels for GeneWeaver using OpenCL."""

import numpy as np
import pyopencl as cl


VECTOR_ADD_KERNEL = """
__kernel void vector_add(
    __global const float *a,
    __global const float *b,
    __global float *result,
    const int n
)
{
    int i = get_global_id(0);

    if (i < n)
    {
        result[i] = a[i] + b[i];
    }
}
"""


def get_opencl_device():
    """Find and return the first available OpenCL GPU."""
    for platform in cl.get_platforms():
        devices = platform.get_devices(
            device_type=cl.device_type.GPU
        )

        if devices:
            return devices[0]

    raise RuntimeError("No OpenCL GPU device available.")


def vector_add_gpu(a, b):
    """
    Add two float32 arrays using an OpenCL GPU kernel.

    Data flow:
        CPU RAM -> GPU memory -> GPU computation -> CPU RAM
    """
    a = np.ascontiguousarray(a, dtype=np.float32)
    b = np.ascontiguousarray(b, dtype=np.float32)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("Input arrays must be one-dimensional.")

    if a.shape != b.shape:
        raise ValueError("Input arrays must have equal lengths.")

    if len(a) == 0:
        return np.array([], dtype=np.float32)

    device = get_opencl_device()

    context = cl.Context(devices=[device])
    queue = cl.CommandQueue(context)

    program = cl.Program(
        context,
        VECTOR_ADD_KERNEL
    ).build()

    result = np.empty_like(a)

    # Allocate GPU memory and copy input arrays to the device.
    a_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=a,
    )

    b_buffer = cl.Buffer(
        context,
        cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=b,
    )

    result_buffer = cl.Buffer(
        context,
        cl.mem_flags.WRITE_ONLY,
        result.nbytes,
    )

    # Launch GPU kernel.
    program.vector_add(
        queue,
        (len(a),),
        None,
        a_buffer,
        b_buffer,
        result_buffer,
        np.int32(len(a)),
    )

    # Copy GPU result back to CPU memory.
    cl.enqueue_copy(
        queue,
        result,
        result_buffer,
    ).wait()

    return result
