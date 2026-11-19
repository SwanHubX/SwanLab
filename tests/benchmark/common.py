"""
@author: nexisato
@file: common.py
@description: bench 公共工具：分位统计、进程查找、内存/CPU 采样、子进程环境、结果落盘。

约束：
  - 只用标准库，不引入 psutil 等依赖；
  - 不修改生产代码，只用现成入口：helper.get_core_impl() 判断实际使用的 core，
    /tmp/swanlab-core-*/core.port 定位 Go 进程，run._core.get_operation_stats() 读上传状态；
  - 路径以仓库根为基准，脚本可在任意 cwd 下运行。
"""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "API_KEY",
    "REPO_ROOT",
    "RESULTS_DIR",
    "BENCH_ROOT",
    "MOCK_SCRIPT",
    "CORE_BINARY",
    "Check",
    "Checks",
    "MockCloud",
    "Sampler",
    "backend_kind",
    "build_child_env",
    "core_pids",
    "environment_info",
    "latency_stats",
    "parse_ps_time",
    "read_port_file_pid",
    "percentile",
    "cleanup_stale_cores",
    "runtime_dir_for_pid",
    "stderr_for_pid",
    "write_result",
    "render_table",
    "scan_local_records",
    "isolation_problems",
    "StatsProbe",
]

BENCH_ROOT = Path(__file__).resolve().parent
REPO_ROOT = BENCH_ROOT.parents[1]
MOCK_SCRIPT = BENCH_ROOT / "mock_cloud" / "server.py"
RESULTS_DIR = BENCH_ROOT / "results"
CORE_BINARY = REPO_ROOT / "swanlab" / "bin" / "swanlab-core"

# Mock Cloud 与测试子进程共用的假凭证；仅在本机回环地址上生效。
API_KEY = "bench-local-api-key"

# 子进程向父进程汇报阶段边界的前缀，父进程据此切分内存采样窗口。
PHASE_MARK = "##PHASE## "
RESULT_MARK = "##BENCH_JSON## "

# --- 子进程环境：阻断所有指向真实后端的配置 ---


def build_child_env(
    mock_url: str,
    root: Optional[Path] = None,
    cwd: Optional[Path] = None,
    core_impl: Optional[str] = None,
    extra: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """构造被测子进程的环境：只允许访问本机 mock，且与宿主凭证完全隔离。

    SWANLAB_ROOT 只覆盖用户级配置目录，cwd 下的 .swanlab/.netrc 需要把 cwd 指到
    临时目录才能隔离，两者都要设置。
    """
    env = {k: v for k, v in os.environ.items() if not k.startswith("SWANLAB_")}
    # 清空所有可能的外联通道
    for key in (
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "ALL_PROXY",
        "http_proxy",
        "https_proxy",
        "all_proxy",
        "SWANLAB_WEBHOOK_URL",
        "SWANLAB_CONFIG_DIR",
        "SWANLAB_SECRETS_DIR",
    ):
        env.pop(key, None)
    env.update(
        {
            "SWANLAB_API_HOST": mock_url,
            "SWANLAB_WEB_HOST": mock_url,
            "SWANLAB_API_KEY": API_KEY,
            "SWANLAB_ON_ONLINE_ERROR": "raise",
            "SWANLAB_ROOT": str(root or (BENCH_ROOT / ".tmp-root")),
            "SWANLAB_CORE_BINARY": str(CORE_BINARY),
            "NO_PROXY": "127.0.0.1,localhost",
            "no_proxy": "127.0.0.1,localhost",
            # 关掉可能拖慢 init 的网络/代理探测
            "PYTHONUNBUFFERED": "1",
        }
    )
    if cwd is not None:
        env["PWD"] = str(cwd)
    if core_impl == "python":
        env["SWANLAB_REQUIRE"] = "core_python"
    elif core_impl == "go":
        env["SWANLAB_REQUIRE"] = "core"
    if extra:
        env.update(extra)
    return env


def isolation_problems() -> List[str]:
    """在子进程内回读解析后的配置，返回所有仍指向真实后端的项（空列表表示隔离成功）。"""
    from swanlab.sdk.internal.settings import Settings

    settings = Settings()
    violations: List[str] = []
    for name, value in (("api_host", settings.api_host), ("web_host", settings.web_host)):
        text = str(value or "")
        if "127.0.0.1" not in text and "localhost" not in text:
            violations.append(f"{name}={text!r} 未指向本机 mock")
    if settings.api_key != API_KEY:
        violations.append(f"api_key 未使用注入值（实际 {settings.api_key!r}）")
    if str(settings.on_online_error) != "raise":
        violations.append(f"on_online_error={settings.on_online_error!r}，失败时不会报错")
    return violations


def backend_kind() -> str:
    """在子进程内调用：返回实际生效的 core 实现（"go" / "python"）。"""
    from swanlab.sdk.internal.pkg import helper

    return helper.get_core_impl()


# --- 统计 ---------------------------------------------------------------------


def percentile(sorted_data: Sequence[float], pct: float) -> float:
    if not sorted_data:
        return 0.0
    idx = int(len(sorted_data) * pct / 100.0)
    return sorted_data[min(idx, len(sorted_data) - 1)]


def latency_stats(values: Sequence[float]) -> Dict[str, float]:
    """微秒延迟的均值与分位；空输入返回全 0。"""
    if not values:
        return {"mean": 0.0, "p50": 0.0, "p95": 0.0, "p99": 0.0}
    ordered = sorted(values)
    return {
        "mean": round(sum(ordered) / len(ordered), 2),
        "p50": round(percentile(ordered, 50), 2),
        "p95": round(percentile(ordered, 95), 2),
        "p99": round(percentile(ordered, 99), 2),
    }


# --- 进程发现 -----------------------------------------------------------------

_PS_TIME_RE = re.compile(r"^(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+(?:\.\d+)?)$")


def parse_ps_time(text: str) -> float:
    """解析 ``ps -o time=`` 的 ``[[dd-]hh:]mm:ss[.cc]`` 为秒。"""
    match = _PS_TIME_RE.match(text.strip())
    if not match:
        return 0.0
    days, hours, minutes, seconds = match.groups()
    return int(days or 0) * 86400 + int(hours or 0) * 3600 + int(minutes or 0) * 60 + float(seconds)


def _run_ps(columns: str) -> List[List[str]]:
    """执行一次 ``ps``，返回按空白切分的行；失败返回空列表。"""
    try:
        out = subprocess.run(
            ["ps", "-axo", columns],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    rows: List[List[str]] = []
    for line in out.stdout.splitlines()[1:]:
        parts = line.split(None, 2)
        if len(parts) >= 2:
            rows.append(parts)
    return rows


def list_processes() -> List[Tuple[int, int, str]]:
    """返回 ``(pid, ppid, command)`` 三元组。"""
    procs: List[Tuple[int, int, str]] = []
    for row in _run_ps("pid=,ppid=,command="):
        try:
            pid, ppid = int(row[0]), int(row[1])
        except ValueError:
            continue
        command = row[2] if len(row) > 2 else ""
        procs.append((pid, ppid, command))
    return procs


def _runtime_dirs() -> List[Path]:
    root = Path("/tmp")
    if not root.is_dir():
        return []
    return sorted(p for p in root.glob("swanlab-core-*") if p.is_dir())


def read_port_file_pid(path: Path) -> Optional[int]:
    """解析 port-file 的 ``pid=`` 行；文件缺失或格式不完整返回 None。"""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if line.startswith("pid="):
            try:
                return int(line[4:])
            except ValueError:
                return None
    return None


def runtime_dir_for_pid(pid: int) -> Optional[Path]:
    """按 core pid 反查其 runtime 目录（内含 port-file 与 stderr）。"""
    for directory in _runtime_dirs():
        if read_port_file_pid(directory / "core.port") == pid:
            return directory
    return None


def stderr_for_pid(pid: int) -> str:
    """读取 core 的 stderr 末尾 8 KiB，用于判断本轮是否出错。"""
    directory = runtime_dir_for_pid(pid)
    if directory is None:
        return ""
    path = directory / "core.stderr.log"
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - 8192))
            return handle.read().decode("utf-8", "replace")
    except OSError:
        return ""


def core_pids(parent_pid: Optional[int] = None) -> List[int]:
    """返回存活的 swanlab-core 进程 pid；``parent_pid`` 非空时只返回其子进程。"""
    found: List[int] = []
    for pid, ppid, command in list_processes():
        if "swanlab-core" not in command:
            continue
        if parent_pid is not None and ppid != parent_pid:
            continue
        found.append(pid)
    return found


def cleanup_stale_cores(timeout: float = 10.0) -> List[int]:
    """结束残留的 core 进程并清理独占模式的临时目录。

    detach=True 的 core 默认要空闲 300s 才退出，残留实例会带着上一轮的指标状态
    和内存继续服务，必须先清掉再开测。
    """
    pids = core_pids()
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not core_pids():
            break
        time.sleep(0.1)
    for pid in core_pids():
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError:
            pass
    for directory in _runtime_dirs():
        if "-own-" in directory.name:
            shutil.rmtree(directory, ignore_errors=True)
    return pids


# --- Mock Cloud ---------------------------------------------------------------


class MockCloud:
    """压测专用 Mock Cloud 子进程。

    独立进程而非线程：既避免与 Python SDK 抢 GIL，也让它不出现在被测进程树的
    内存统计里（它的开销单独采样上报）。
    """

    def __init__(self, workdir: Path, latency_ms: int = 0) -> None:
        self.workdir = workdir
        self.latency_ms = latency_ms
        self.port_file = workdir / "mock.port"
        self.log_path = workdir / "mock.log"
        self.process: Optional[subprocess.Popen] = None
        self.port = 0
        self._log: Optional[Any] = None

    def __enter__(self) -> "MockCloud":
        self.workdir.mkdir(parents=True, exist_ok=True)
        self._log = self.log_path.open("w", encoding="utf-8")
        argv = [sys.executable, str(MOCK_SCRIPT), "--port", "0", "--port-file", str(self.port_file)]
        if self.latency_ms:
            argv.extend(["--latency-ms", str(self.latency_ms)])
        self.process = subprocess.Popen(
            argv,
            cwd=str(self.workdir),
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
            stdout=self._log,
            stderr=subprocess.STDOUT,
        )
        self.port = self._wait_port()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=5)
        if self._log is not None:
            self._log.close()

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.port}"

    @property
    def pid(self) -> Optional[int]:
        return self.process.pid if self.process is not None else None

    def _wait_port(self, timeout: float = 20.0) -> int:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.process is not None and self.process.poll() is not None:
                raise RuntimeError(f"Mock Cloud 提前退出，日志见 {self.log_path}")
            if self.port_file.exists():
                text = self.port_file.read_text(encoding="utf-8").strip()
                if text:
                    return int(text)
            time.sleep(0.02)
        raise RuntimeError(f"Mock Cloud 端口文件未就绪：{self.port_file}")

    def stats(self) -> Dict[str, Any]:
        """拉取 mock 的计数快照（只含计数与体积，不含任何记录明细）。"""
        import urllib.request

        with urllib.request.urlopen(f"{self.url}/__mock/stats", timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))


# --- 内存 / CPU 采样 -----------------------------------------------------------

Sample = Dict[str, Any]


class Sampler(threading.Thread):
    """按固定间隔采集一组 pid 的常驻内存与累计 CPU 时间。

    CPU 取 ps -o time= 的累计秒数差，macOS 上分辨率是 1 秒；内存取 ps 的 RSS
    （不含被换出或共享计费差异，仅用于横向对比）。
    """

    def __init__(self, pids_fn: Callable[[], Iterable[int]], interval: float = 0.1) -> None:
        super().__init__(daemon=True, name="bench-sampler")
        self._pids_fn = pids_fn
        self._interval = interval
        self._halt = threading.Event()
        self.samples: List[Sample] = []

    def run(self) -> None:
        while not self._halt.is_set():
            pids = sorted({int(p) for p in self._pids_fn() if p})
            self.samples.append(self._sample(pids))
            self._halt.wait(self._interval)

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=5)

    @staticmethod
    def _sample(pids: List[int]) -> Sample:
        started = time.monotonic()
        rss: Dict[int, int] = {}
        cpu: Dict[int, float] = {}
        if pids:
            try:
                out = subprocess.run(
                    ["ps", "-o", "pid=,rss=,time=", "-p", ",".join(str(p) for p in pids)],
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                for line in out.stdout.splitlines():
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        pid = int(parts[0])
                    except ValueError:
                        continue
                    rss[pid] = int(parts[1]) * 1024
                    cpu[pid] = parse_ps_time(parts[2])
            except (OSError, subprocess.SubprocessError):
                pass
        return {"t": started, "rss": rss, "cpu": cpu}

    # -- 查询辅助 --------------------------------------------------------------

    def peak_rss(self, pids: Sequence[int], start: float = 0.0, end: float = 1e18) -> int:
        peak = 0
        for sample in self.samples:
            if not (start <= sample["t"] <= end):
                continue
            peak = max(peak, sum(sample["rss"].get(pid, 0) for pid in pids))
        return peak

    def cpu_seconds(self, pids: Sequence[int]) -> float:
        """各 pid 的累计 CPU 时间之和；进程退出后不再更新。"""
        totals: Dict[int, float] = {}
        for sample in self.samples:
            for pid in pids:
                value = sample["cpu"].get(pid)
                if value is not None:
                    totals[pid] = max(totals.get(pid, 0.0), value)
        return sum(totals.values())


# --- 校验与报告 ---------------------------------------------------------------


class Check:
    """一条校验结果。"""

    def __init__(self, name: str, ok: bool, detail: str = "") -> None:
        self.name = name
        self.ok = ok
        self.detail = detail

    def __repr__(self) -> str:  # pragma: no cover - 调试输出
        return f"Check({self.name!r}, ok={self.ok}, detail={self.detail!r})"


class Checks:
    def __init__(self) -> None:
        self.items: List[Check] = []

    def add(self, name: str, ok: bool, detail: str = "") -> None:
        self.items.append(Check(name, ok, detail))

    def equal(self, name: str, actual: Any, expected: Any) -> None:
        self.add(name, actual == expected, f"期望 {expected!r}，实际 {actual!r}")

    def ok(self) -> bool:
        return all(item.ok for item in self.items)

    def failures(self) -> List[Check]:
        return [item for item in self.items if not item.ok]


def environment_info() -> Dict[str, Any]:
    """收集解释跨天、跨机差异所需的环境信息。"""
    import platform

    def _git(*args: str) -> str:
        try:
            return subprocess.run(
                ["git", *args], cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=15
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""

    binary_sha = ""
    if CORE_BINARY.exists():
        import hashlib

        binary_sha = hashlib.sha256(CORE_BINARY.read_bytes()).hexdigest()

    return {
        "git_sha": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "core_binary": str(CORE_BINARY),
        "core_binary_sha256": binary_sha,
        "machine": platform.machine(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu_count": os.cpu_count(),
        "GOGC": os.environ.get("GOGC", ""),
        "GOMEMLIMIT": os.environ.get("GOMEMLIMIT", ""),
    }


def write_result(payload: Dict[str, Any], name: str) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    latest = RESULTS_DIR / "latest.json"
    latest.write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    return path


def render_table(headers: Sequence[str], rows: Sequence[Sequence[Any]], widths: Optional[Sequence[int]] = None) -> str:
    """左对齐制表；首列为标签列。"""
    cols = widths or [max(len(str(headers[i])), *(len(str(r[i])) for r in rows)) + 2 for i in range(len(headers))]
    lines = ["  " + "".join(str(headers[i]).ljust(cols[i]) for i in range(len(headers)))]
    lines.append("  " + "-" * (sum(cols) + 2))
    for row in rows:
        lines.append("  " + "".join(str(row[i]).ljust(cols[i]) for i in range(len(headers))))
    return "\n".join(lines)


# --- 子进程内的采集与回读 -----------------------------------------------------


def scan_local_records(run_dir: Path) -> Dict[str, int]:
    """回读 run-*.swanlab，按 record 类型计数。"""
    from swanlab.proto.swanlab.record.v1 import record_pb2
    from swanlab.sdk.internal.core_python.store import DataStoreReader

    counts: Dict[str, int] = {}
    for path in sorted(Path(run_dir).glob("run-*.swanlab")):
        reader = DataStoreReader()
        reader.open(path)
        try:
            for raw in reader:
                record = record_pb2.Record()
                record.ParseFromString(raw)
                kind = record.WhichOneof("record_type") or "unknown"
                counts[kind] = counts.get(kind, 0) + 1
        finally:
            reader.close()
    return counts


class StatsProbe(threading.Thread):
    """周期性读取 core 的 OperationStats（在子进程内运行，因为需要 Run 对象）。

    间隔 0.2s 而不是 1s：大场景运行时间只有几十秒，1s 会漏掉 finish 阶段的峰值，
    而单次调用只是一次本地 gRPC。
    """

    def __init__(self, run: Any, interval: float = 0.2) -> None:
        super().__init__(daemon=True, name="bench-stats")
        self._run = run
        self._interval = interval
        self._halt = threading.Event()
        self.samples: List[Dict[str, Any]] = []
        self.errors: List[str] = []

    def run(self) -> None:
        while not self._halt.is_set():
            self._poll()
            self._halt.wait(self._interval)

    def _poll(self) -> None:
        try:
            response = self._run._core.get_operation_stats()
        except Exception as exc:  # noqa: BLE001 - 采样失败不应影响被测链路
            self.errors.append(f"{type(exc).__name__}: {exc}")
            return
        if not response.success or not response.HasField("stats"):
            return
        stats = response.stats
        upload = stats.upload
        self.samples.append(
            {
                "paused": bool(upload.paused),
                "buffered_bytes": int(upload.buffered_bytes),
                "pending_records": int(upload.pending_records),
                "uploaded_records": int(stats.uploaded_records),
                "state": int(stats.state),
            }
        )

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=5)

    def summary(self) -> Dict[str, Any]:
        return {
            "paused_seen": any(sample["paused"] for sample in self.samples),
            "max_buffered_bytes": max((s["buffered_bytes"] for s in self.samples), default=0),
            "max_pending_records": max((s["pending_records"] for s in self.samples), default=0),
            "max_uploaded_records": max((s["uploaded_records"] for s in self.samples), default=0),
            "samples": len(self.samples),
            "errors": self.errors[:5],
        }
