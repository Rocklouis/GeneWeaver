
from types import SimpleNamespace

import pyopencl as cl

from geneweaver.gpu.device import (
    detect_gpu,
    get_compute_device,
    print_device_info,
)


def test_detect_gpu():
    fake_device = SimpleNamespace(
        name="Test GPU",
        global_mem_size=1024 * 1024 * 512,
        max_work_group_size=256,
    )

    fake_platform = SimpleNamespace(
        get_devices=lambda device_type: [fake_device]
    )

    original_get_platforms = cl.get_platforms
    cl.get_platforms = lambda: [fake_platform]

    try:
        result = detect_gpu()

        assert result.backend == "OpenCL"
        assert result.name == "Test GPU"
        assert result.is_available is True
        assert result.device_type == "GPU"
        assert result.global_memory_mb == 512
        assert result.max_work_group_size == 256
    finally:
        cl.get_platforms = original_get_platforms


def test_cpu_fallback_when_gpu_unavailable(monkeypatch):
    monkeypatch.setattr(
        cl,
        "get_platforms",
        lambda: (_ for _ in ()).throw(RuntimeError("GPU unavailable")),
    )

    result = detect_gpu()

    assert result.backend == "CPU"
    assert result.name == "CPU Fallback"
    assert result.is_available is False
    assert result.device_type == "CPU"


def test_cpu_fallback_when_no_gpu_found(monkeypatch):
    fake_platform = SimpleNamespace(
        get_devices=lambda device_type: []
    )

    monkeypatch.setattr(
        cl,
        "get_platforms",
        lambda: [fake_platform],
    )

    result = detect_gpu()

    assert result.backend == "CPU"
    assert result.is_available is False


def test_get_compute_device():
    result = get_compute_device()

    assert result.backend in ("OpenCL", "CPU")
    assert result.device_type in ("GPU", "CPU")


def test_print_device_info(capsys):
    result = detect_gpu()
    print_device_info(result)

    output = capsys.readouterr().out

    assert "GENEWEAVER DEVICE INFORMATION" in output
    assert f"Backend: {result.backend}" in output
    assert f"Device: {result.name}" in output
