"""
@author: cunyue
@description: CLI ping 模块，连通性与环境诊断命令

预期输出样式：

╭─ SwanLab Diagnostics ─────────────────────────────────────╮
│ Server link: OK · https://api.swanlab.cn · 42 ms           │
│                                                           │
│ 1. SDK                                                    │
│   ✓ Version        0.7.0                                  │
│   ✓ Python         3.11.8                                 │
│   ✓ Mode           online                                 │
│                                                           │
│ 2. System                                                 │
│   ✓ OS             macOS 15.5 arm64                       │
│   ✓ Python         3.11.8                                 │
│   ✓ Executable     /path/to/python                        │
│                                                           │
│ 3. Hardware                                               │
│   ✓ Apple Silicon  Apple M4 Pro                           │
│   ✓ CPU Cores      14                                     │
│   ✓ Unified Memory 48 GB                                  │
│                                                           │
│   Accelerators                                            │
│     - Not detected                                        │
╰───────────────────────────────────────────────────────────╯

硬件展示说明：
  - Accelerators 不要展示成单个 GPU 字段。probe_python 中加速器是 list[AcceleratorSnapshot]，
    需要支持多卡、多厂商异构场景。
  - 有加速器时按 vendor 分组，再列出每个 device，单卡也沿用同一结构，避免输出形状随设备数量变化。

异构加速器示例：
  Accelerators
    ✓ NVIDIA        2 devices, driver 550.54, CUDA 12.4
      [0] RTX 4090  24 GB
      [1] RTX 4090  24 GB
    ✓ Ascend        8 devices, CANN 8.0
      [0] Ascend 910B 64 GB

状态标识：
  ✓ 正常
  ! 警告
  ✗ 异常
  - 未检测到 / 未配置 / 不适用
"""

import asyncio
import sys
import time
from dataclasses import dataclass
from typing import Optional, Union

import click
import requests
from google.protobuf.message import Message
from rich.console import Console
from rich.markup import escape
from rich.panel import Panel

from swanlab.proto.swanlab.env.v1.metadata_pb2 import (
    AppleSiliconSnapshot,
    DeviceSnapshot,
    HardwareSnapshot,
    MemorySnapshot,
    MetadataSnapshot,
)
from swanlab.proto.swanlab.grpc.probe.v1.probe_pb2 import DeliverProbeStartRequest, GetMetadataSnapshotResponse
from swanlab.sdk import Settings, impl, pkg
from swanlab.sdk.internal.pkg import adapter

MemoryLikeSnapshot = Union[AppleSiliconSnapshot, DeviceSnapshot, MemorySnapshot]


@click.command()
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host of the swanlab server.",
)
def ping(host: Optional[str]):
    """Diagnose connectivity and system environment for the current machine."""
    args = {"api_host": host}
    settings = Settings(**pkg.helper.strip_none(args))

    # 定义变量来存储检查结果，供后续展示使用
    check_server_result: Optional[CheckServerResult] = None
    check_probe_result: Optional[GetMetadataSnapshotResponse] = None

    # 1. 调用接口，检查服务器连通性和响应状态
    async def check_server():
        nonlocal check_server_result

        endpoint: str = settings.api_host
        api_host = endpoint + "/api"

        def request_server() -> CheckServerResult:
            start = time.perf_counter()
            try:
                session = pkg.client.session.create(timeout=5, default_retry=0)
                response = session.get(api_host)
                latency_ms = int((time.perf_counter() - start) * 1000)
                if response.ok:
                    text = response.text.strip()
                    if text != "OK":
                        return CheckServerResult(
                            ok=False,
                            status=f"UNEXPECTED RESPONSE: {text}",
                            endpoint=endpoint,
                            latency_ms=latency_ms,
                        )
                    return CheckServerResult(ok=True, status="OK", endpoint=endpoint, latency_ms=latency_ms)

                return CheckServerResult(
                    ok=False,
                    status=f"HTTP {response.status_code}",
                    endpoint=endpoint,
                    latency_ms=latency_ms,
                )
            except requests.Timeout:
                return CheckServerResult(ok=False, status="TIMEOUT", endpoint=endpoint, latency_ms=5000)
            except requests.RequestException:
                latency_ms = int((time.perf_counter() - start) * 1000)
                return CheckServerResult(ok=False, status="ERROR", endpoint=endpoint, latency_ms=latency_ms)

        check_server_result = await asyncio.to_thread(request_server)

    # 2. 启动probe获取当前环境信息
    async def check_probe():
        nonlocal check_probe_result

        # probe 被设计为，如果 core 不可访问，依旧采集信息但不上报，我们可以通过对应的rpc函数拿到 probe 采集的环境信息
        probe = impl.create_probe()
        probe.deliver_probe_start(DeliverProbeStartRequest(probe_settings=settings.to_probe_proto()))
        check_probe_result = probe.get_metadata_snapshot()

    # 3. 并行执行服务器检查和probe检查，减少等待时间
    async def main():
        await asyncio.gather(
            check_server(),
            check_probe(),
        )

    asyncio.run(main())

    show_result(settings, check_server_result, check_probe_result)


@dataclass
class CheckServerResult:
    ok: bool
    status: str
    endpoint: str
    latency_ms: int


def show_result(
    settings: Settings,
    server_result: Optional[CheckServerResult],
    probe_result: Optional[GetMetadataSnapshotResponse],
):
    """展示检查结果，格式化输出"""
    console = Console()

    # 1. 校验检查结果，失败时给出明确的排查入口
    if server_result is None:
        console.print("[red]Failed to check server connectivity.[/red] Use `SWANLAB_DEBUG=1 swanlab ping` for details.")
        sys.exit(1)
    if probe_result is None:
        console.print("[red]Failed to collect probe information.[/red] Use `SWANLAB_DEBUG=1 swanlab ping` for details.")
        sys.exit(1)
    if not probe_result.success:
        console.print(
            f"[red]Probe reported failure:[/red] {escape(probe_result.message)}. "
            "Use `SWANLAB_DEBUG=1 swanlab ping` for details."
        )
        sys.exit(1)
    if not probe_result.HasField("metadata"):
        console.print("[red]Probe did not return metadata.[/red] Use `SWANLAB_DEBUG=1 swanlab ping` for details.")
        sys.exit(1)

    metadata = probe_result.metadata
    status_style = "green" if server_result.ok else "red"

    # 2. 按用户可读的诊断分组组装展示内容
    server_link = f"Server {server_result.status} · {server_result.endpoint} · {server_result.latency_ms} ms"
    lines = [
        f"[{status_style}]{escape(server_link)}[/{status_style}]",
        "",
        "[bold cyan]1. SDK[/bold cyan]",
        _result_line(True, "Version", pkg.helper.get_swanlab_version()),
        _result_line(True, "Python", _metadata_python_version(metadata) or sys.version.split()[0]),
        _result_line(True, "Mode", settings.mode),
        "",
        "[bold cyan]2. System[/bold cyan]",
    ]
    lines.extend(_format_system(metadata))
    lines.append("")
    lines.append("[bold cyan]3. Hardware[/bold cyan]")
    lines.extend(_format_hardware(metadata))

    # 3. 使用 rich 输出最终 UI
    console.print(Panel("\n".join(lines), title="SwanLab Diagnostics", border_style=status_style))


def _has_field(message: Message, field: str) -> bool:
    try:
        return message.HasField(field)
    except ValueError:
        return False


def _mark(ok: bool) -> str:
    return "[green]✓[/green]" if ok else "[red]✗[/red]"


def _missing() -> str:
    return "[dim]-[/dim]"


def _result_line(ok: bool, label: str, value: str) -> str:
    return f"  {_mark(ok)} {label:<14} {escape(value)}"


def _missing_line(label: str, value: str = "Not detected") -> str:
    return f"  {_missing()} {label:<14} [dim]{escape(value)}[/dim]"


def _optional_string(message: Message, field: str) -> Optional[str]:
    if not _has_field(message, field):
        return None
    value = getattr(message, field)
    if value == "":
        return None
    return str(value)


def _metadata_python_version(metadata: MetadataSnapshot) -> Optional[str]:
    if not _has_field(metadata, "runtime"):
        return None
    return _optional_string(metadata.runtime, "python_version")


def _format_memory(message: MemoryLikeSnapshot, value_field: str, unit_field: str) -> Optional[str]:
    if not _has_field(message, value_field):
        return None
    value = getattr(message, value_field)
    if _has_field(message, unit_field):
        unit = adapter.memory_unit.get(getattr(message, unit_field))
        if isinstance(unit, str):
            return f"{value} {unit}"
    return str(value)


def _format_vendor(vendor: int) -> str:
    value = adapter.accelerator_vendor.get(vendor)
    if not isinstance(value, str):
        return "Unknown"
    if value == "rocm":
        return "ROCm"
    if value == "nvidia":
        return "NVIDIA"
    return value.title()


def _format_system(metadata: MetadataSnapshot) -> list[str]:
    if not _has_field(metadata, "runtime"):
        return [_missing_line("Runtime", "Not collected")]

    runtime = metadata.runtime
    lines = []
    os_name = _optional_string(runtime, "os_pretty") or _optional_string(runtime, "os")
    lines.append(_result_line(True, "OS", os_name) if os_name else _missing_line("OS", "Not collected"))
    python_version = _optional_string(runtime, "python_version")
    lines.append(
        _result_line(True, "Python", python_version) if python_version else _missing_line("Python", "Not collected")
    )
    executable = _optional_string(runtime, "python_executable")
    lines.append(
        _result_line(True, "Executable", executable) if executable else _missing_line("Executable", "Not collected")
    )
    return lines


def _format_hardware(metadata: MetadataSnapshot) -> list[str]:
    if not _has_field(metadata, "hardware"):
        return [_missing_line("Hardware", "Not collected")]

    hardware = metadata.hardware
    lines = []
    if _has_field(hardware, "apple_silicon"):
        apple = hardware.apple_silicon
        name = _optional_string(apple, "name")
        lines.append(
            _result_line(True, "Apple Silicon", name) if name else _missing_line("Apple Silicon", "Not detected")
        )
        cpu_count = _optional_string(apple, "cpu_count")
        lines.append(_result_line(True, "CPU Cores", cpu_count) if cpu_count else _missing_line("CPU Cores"))
        memory = _format_memory(apple, "memory", "memory_unit")
        lines.append(_result_line(True, "Unified Memory", memory) if memory else _missing_line("Unified Memory"))
    else:
        lines.extend(_format_cpu_memory(hardware))

    lines.append("")
    lines.append("  [bold]Accelerators[/bold]")
    lines.extend(_format_accelerators(hardware))
    return lines


def _format_cpu_memory(hardware: HardwareSnapshot) -> list[str]:
    lines = []
    if _has_field(hardware, "cpu"):
        cpu = hardware.cpu
        brand = _optional_string(cpu, "brand")
        lines.append(_result_line(True, "CPU", brand) if brand else _missing_line("CPU"))
        core_parts = []
        if _has_field(cpu, "logical_count"):
            core_parts.append(f"{cpu.logical_count} logical")
        if _has_field(cpu, "physical_count"):
            core_parts.append(f"{cpu.physical_count} physical")
        lines.append(
            _result_line(True, "CPU Cores", " / ".join(core_parts)) if core_parts else _missing_line("CPU Cores")
        )
    else:
        lines.append(_missing_line("CPU", "Not collected"))
        lines.append(_missing_line("CPU Cores", "Not collected"))

    if _has_field(hardware, "memory"):
        memory = _format_memory(hardware.memory, "total", "total_unit")
        lines.append(_result_line(True, "Memory", memory) if memory else _missing_line("Memory"))
    else:
        lines.append(_missing_line("Memory", "Not collected"))
    return lines


def _format_accelerators(hardware: HardwareSnapshot) -> list[str]:
    if not hardware.accelerators:
        return [f"    {_missing()} [dim]Not detected[/dim]"]

    lines = []
    for accelerator in hardware.accelerators:
        vendor = _format_vendor(accelerator.vendor)
        details = [f"{len(accelerator.devices)} devices"]
        if _has_field(accelerator, "version"):
            details.append(f"driver {accelerator.version}")
        if _has_field(accelerator, "cuda_version"):
            details.append(f"CUDA {accelerator.cuda_version}")
        if _has_field(accelerator, "cann_version"):
            details.append(f"CANN {accelerator.cann_version}")
        lines.append(f"    {_mark(True)} {vendor:<14} {escape(', '.join(details))}")
        for device in accelerator.devices:
            index = str(device.index) if _has_field(device, "index") else "-"
            name = _optional_string(device, "name") or "Unknown"
            memory = _format_memory(device, "memory", "memory_unit")
            suffix = f"  {memory}" if memory else ""
            lines.append(f"      {escape(f'[{index}] {name}{suffix}')}")
    return lines
