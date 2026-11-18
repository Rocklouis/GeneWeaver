
"""GPU device detection and CPU fallback for GeneWeaver."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class DeviceInfo:
    backend: str
    name: str
    is_available: bool
    device_type: str
    global_memory_mb: Optional[float] = None
    max_work_group_size: Optional[int] = None


def detect_gpu() -> DeviceInfo:
    """
    Detect an available GPU using OpenCL.

    Returns GPU information when an OpenCL GPU is available.
    Otherwise, returns CPU fallback information.
    """
    try:
        import pyopencl as cl

        platforms = cl.get_platforms()

        for platform in platforms:
            for device in platform.get_devices(
                device_type=cl.device_type.GPU
            ):
                return DeviceInfo(
                    backend="OpenCL",
                    name=device.name.strip(),
                    is_available=True,
                    device_type="GPU",
                    global_memory_mb=round(
                        device.global_mem_size / (1024 * 1024), 2
                    ),
                    max_work_group_size=device.max_work_group_size,
                )

    except Exception as exc:
        print(f"GPU detection warning: {exc}")

    return DeviceInfo(
        backend="CPU",
        name="CPU Fallback",
        is_available=False,
        device_type="CPU",
    )


def print_device_info(info: DeviceInfo) -> None:
    """Display detected processing device information."""
    print("=" * 45)
    print("GENEWEAVER DEVICE INFORMATION")
    print("=" * 45)
    print(f"Backend: {info.backend}")
    print(f"Device: {info.name}")
    print(f"Available: {info.is_available}")
    print(f"Device type: {info.device_type}")

    if info.is_available:
        print(f"Global memory: {info.global_memory_mb} MB")
        print(f"Maximum work group size: {info.max_work_group_size}")

    print("=" * 45)


def get_compute_device() -> DeviceInfo:
    """Return the detected GPU or CPU fallback."""
    return detect_gpu()
