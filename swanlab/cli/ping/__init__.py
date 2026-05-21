"""
@author: cunyue
@description: CLI ping 模块，连通性与环境诊断命令

预期输出样式：

╭─ SwanLab Diagnostics ─────────────────────────────────────╮
│ Status: OK                                                │
│ Endpoint: https://api.swanlab.cn                          │
╰───────────────────────────────────────────────────────────╯

SDK
  ✓ Version        0.7.0
  ✓ Python         3.11.8
  ✓ Mode           cloud

Server
  ✓ Endpoint       https://api.swanlab.cn
  ✓ Status         OK
  ✓ Latency        42 ms

System
  ✓ OS             macOS 15.5 arm64
  ✓ Python         3.11.8
  ✓ Executable     /path/to/python

Hardware
  ✓ Apple Silicon  Apple M4 Pro
  ✓ CPU Cores      14
  ✓ Unified Memory 48 GB

  Accelerators
    - Not detected

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
from typing import Optional

import click
import requests

from swanlab.proto.swanlab.grpc.probe.v1.probe_pb2 import DeliverProbeStartRequest, GetMetadataSnapshotResponse
from swanlab.sdk import Settings, impl, pkg


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

    show_result(check_server_result, check_probe_result)


@dataclass
class CheckServerResult:
    ok: bool
    status: str
    endpoint: str
    latency_ms: int


def show_result(server_result: Optional[CheckServerResult], probe_result: Optional[GetMetadataSnapshotResponse]):
    """展示检查结果，格式化输出"""
    if server_result is None:
        click.echo("Failed to check server connectivity, use `SWANLAB_DEBUG=1 swanlab ping` for more details.")
        sys.exit(1)
    if probe_result is None:
        click.echo("Failed to collect probe information, use `SWANLAB_DEBUG=1 swanlab ping` for more details.")
        sys.exit(1)
    if not probe_result.success:
        click.echo(
            f"Probe reported failure with message: {probe_result.message}, use `SWANLAB_DEBUG=1 swanlab ping` for more details."
        )
        sys.exit(1)
    print(server_result)
