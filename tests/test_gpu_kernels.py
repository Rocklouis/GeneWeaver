
import numpy as np
import pytest

from geneweaver.gpu.kernels import vector_add_gpu


def test_vector_add_small_arrays():
    a = np.array([1, 2, 3, 4], dtype=np.float32)
    b = np.array([5, 6, 7, 8], dtype=np.float32)

    result = vector_add_gpu(a, b)

    np.testing.assert_array_equal(result, a + b)


def test_vector_add_large_arrays():
    a = np.arange(100_000, dtype=np.float32)
    b = np.ones(100_000, dtype=np.float32)

    result = vector_add_gpu(a, b)

    np.testing.assert_array_equal(result, a + b)


def test_vector_add_empty_arrays():
    a = np.array([], dtype=np.float32)
    b = np.array([], dtype=np.float32)

    result = vector_add_gpu(a, b)

    assert result.size == 0
    assert result.dtype == np.float32


def test_vector_add_unequal_lengths():
    a = np.array([1, 2, 3], dtype=np.float32)
    b = np.array([4, 5], dtype=np.float32)

    with pytest.raises(ValueError, match="equal lengths"):
        vector_add_gpu(a, b)


def test_vector_add_multidimensional_arrays():
    a = np.array([[1, 2], [3, 4]], dtype=np.float32)
    b = np.array([[5, 6], [7, 8]], dtype=np.float32)

    with pytest.raises(ValueError, match="one-dimensional"):
        vector_add_gpu(a, b)


def test_vector_add_non_contiguous_arrays():
    a = np.arange(10, dtype=np.float32)[::2]
    b = np.ones(5, dtype=np.float32)

    assert not a.flags["C_CONTIGUOUS"]

    result = vector_add_gpu(a, b)

    np.testing.assert_array_equal(result, a + b)


def test_vector_add_device_unavailable(monkeypatch):
    def unavailable_device():
        raise RuntimeError("No OpenCL GPU device available.")

    monkeypatch.setattr(
        "geneweaver.gpu.kernels.get_opencl_device",
        unavailable_device,
    )

    a = np.array([1, 2], dtype=np.float32)
    b = np.array([3, 4], dtype=np.float32)

    with pytest.raises(RuntimeError, match="No OpenCL GPU"):
        vector_add_gpu(a, b)
