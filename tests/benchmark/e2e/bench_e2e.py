"""
@author: nexisato
@file: bench_e2e.py
@description: 端到端基准：online 模式 + 本地 Mock Cloud，对比 No-SDK / Python core / Go core
              在场景 A/B/C 下的初始化耗时、run.log 延迟、吞吐、内存与日志体积。

父进程负责启动 Mock Cloud 并采样内存/CPU；被测代码在子进程里跑，避免全局单例和上一轮
状态互相影响。子进程用 ##PHASE## 行上报阶段边界，父进程据此切分采样窗口。

用法::

    # 冒烟（1 万条）
    uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --steps 100 --keys 100 --core both

    # 正式场景 A/C
    uv run python tests/benchmark/e2e/bench_e2e.py --scenario A --core both --repeats 3
    uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --skip-store
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import common as H  # noqa: E402

# 场景负载；C 为默认规模，A/B 沿用历史报告的规模。
SCENARIOS: Dict[str, Dict[str, int]] = {
    "A": {"steps": 100_000, "keys": 9},
    "B": {"steps": 7, "keys": 100_000},
    "C": {"steps": 1_000, "keys": 100},
}

# 固定的 core 参数。record_batch 同时决定 consumer 小批(<=100)、gRPC 批(<=2000)、
# 上传批大小和 Go 侧队列容量(8x)，改它会同时影响多个环节。
PROFILES: Dict[str, Dict[str, float]] = {
    "default": {"record_batch": 10_000, "record_interval": 1.5},
    "isolated": {"record_batch": 10_000, "record_interval": 3_600.0},
}

CORES = ("none", "python", "go")


# ===========================================================================
# 子进程 Worker
# ===========================================================================


def _emit_phase(name: str, **extra: Any) -> None:
    payload = {"name": name, "t": time.monotonic()}
    payload.update(extra)
    print(H.PHASE_MARK + json.dumps(payload), flush=True)


def run_worker(args: argparse.Namespace) -> None:
    """被测子进程：一次完整的 init → log → finish。"""
    import swanlab

    core_impl = H.backend_kind()
    violations = H.isolation_problems()
    _emit_phase("ready", core_impl=core_impl)

    settings = swanlab.Settings(
        core=swanlab.Settings.Core(
            record_batch=int(args.record_batch),
            record_interval=float(args.record_interval),
            skip_store=bool(args.skip_store),
        ),
        probe=swanlab.Settings.Probe(monitor=False),
    )

    t_init_start = time.perf_counter()
    run = swanlab.init(
        mode="online",
        project=args.project,
        experiment_name=args.experiment,
        settings=settings,
    )
    init_ms = (time.perf_counter() - t_init_start) * 1000.0
    core_pids = H.core_pids(parent_pid=os.getpid())
    _emit_phase("init_done", init_ms=round(init_ms, 2), core_pids=core_pids)

    probe = H.StatsProbe(run)
    probe.start()

    keys = list(range(int(args.keys)))
    latencies: List[float] = []
    producer_s = 0.0
    finish_ms = 0.0
    total_s = 0.0
    local_counts: Dict[str, int] = {}
    run_file_bytes = 0
    error: Optional[str] = None
    media_logged = 0
    save_files = 0
    media_every = max(0, int(args.media_every))
    save_count = max(0, int(args.save_count))
    frame: Any = None
    if media_every > 0:
        import numpy as np

        frame = np.zeros((16, 16, 3), dtype=np.uint8)

    latency_head_us = round(latencies[0], 2) if latencies else 0.0
    run_dir_str = ""
    try:
        # 混合档：先登记待上传文件，policy=end 表示到 finish 才上传
        for index in range(save_count):
            ckpt = Path.cwd() / f"ckpt-{index}.bin"
            ckpt.write_bytes(bytes([index % 256]) * 64 * 1024)
            save_files += len(run.save(str(ckpt), policy="end"))

        t_total_start = time.perf_counter()
        for step in range(int(args.steps)):
            data = {f"metric_{index}": 0.5 + (step % 1000) * 1e-4 for index in keys}
            t0 = time.perf_counter()
            run.log(data, step=step)
            latencies.append((time.perf_counter() - t0) * 1e6)
            if media_every > 0 and step % media_every == 0:
                # 每张图内容不同，文件名含 step 与内容哈希，不会互相覆盖
                frame.fill(step % 256)
                run.log_image(key="image", data=frame, step=step)
                media_logged += 1
            if step % 2000 == 0:
                _emit_phase("progress", step=step)

        t_producer_end = time.perf_counter()
        producer_s = t_producer_end - t_total_start
        _emit_phase("producer_done", producer_s=round(producer_s, 4))

        t_finish_start = time.perf_counter()
        run.finish()
        finish_ms = (time.perf_counter() - t_finish_start) * 1000.0
        t_total_end = time.perf_counter()
        total_s = t_total_end - t_total_start
        _emit_phase("finish_done", finish_ms=round(finish_ms, 2))

        probe.stop()
        run_dir = Path(run._ctx.run_dir)
        run_dir_str = str(run_dir)
        for path in run_dir.glob("run-*.swanlab"):
            run_file_bytes += path.stat().st_size
    except BaseException as exc:  # noqa: BLE001 - 任何异常都要带上下文回传给父进程
        error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        probe.stop()

    result = {
        "ok": error is None,
        "error": error,
        "core_impl": core_impl,
        "isolation_problems": violations,
        "steps": int(args.steps),
        "keys": int(args.keys),
        "skip_store": bool(args.skip_store),
        "record_batch": int(args.record_batch),
        "record_interval": float(args.record_interval),
        "init_ms": round(init_ms, 2),
        "producer_s": round(producer_s, 4),
        "finish_ms": round(finish_ms, 2),
        "total_s": round(total_s, 4),
        "scalar_expected": int(args.steps) * int(args.keys),
        "latency_us": H.latency_stats(latencies),
        "latency_head_us": round(latencies[0], 2) if latencies else latency_head_us,
        "run_dir": run_dir_str,
        "local_records": local_counts,
        "run_file_bytes": run_file_bytes,
        "media_logged": media_logged,
        "save_files": save_files,
        "upload": probe.summary(),
    }
    print(H.RESULT_MARK + json.dumps(result), flush=True)


def _consume(data: Dict[str, float]) -> int:
    return len(data)


def run_worker_none(args: argparse.Namespace) -> None:
    """No-SDK 基线：与 run_worker 同构的循环，但不导入 swanlab。"""
    _emit_phase("ready", core_impl="none")
    t_init_start = time.perf_counter()
    init_ms = (time.perf_counter() - t_init_start) * 1000.0
    core_pids = H.core_pids(parent_pid=os.getpid())
    _emit_phase("init_done", init_ms=round(init_ms, 2), core_pids=core_pids)

    keys = list(range(int(args.keys)))
    latencies: List[float] = []
    t_total_start = time.perf_counter()
    for step in range(int(args.steps)):
        data = {f"metric_{index}": 0.5 + (step % 1000) * 1e-4 for index in keys}
        t0 = time.perf_counter()
        _consume(data)
        latencies.append((time.perf_counter() - t0) * 1e6)
        if step % 2000 == 0:
            _emit_phase("progress", step=step)
    producer_s = time.perf_counter() - t_total_start
    _emit_phase("producer_done", producer_s=round(producer_s, 4))
    total_s = time.perf_counter() - t_total_start
    _emit_phase("finish_done", finish_ms=0.0)

    result = {
        "ok": True,
        "error": None,
        "core_impl": "none",
        "isolation_problems": [],
        "steps": int(args.steps),
        "keys": int(args.keys),
        "skip_store": False,
        "record_batch": int(args.record_batch),
        "record_interval": float(args.record_interval),
        "init_ms": round(init_ms, 2),
        "producer_s": round(producer_s, 4),
        "finish_ms": 0.0,
        "total_s": round(total_s, 4),
        "scalar_expected": int(args.steps) * int(args.keys),
        "latency_us": H.latency_stats(latencies),
        "latency_head_us": round(latencies[0], 2) if latencies else 0.0,
        "local_records": {},
        "run_file_bytes": 0,
        "media_logged": 0,
        "save_files": 0,
        "upload": {
            "paused_seen": False,
            "max_buffered_bytes": 0,
            "max_pending_records": 0,
            "max_uploaded_records": 0,
            "samples": 0,
            "errors": [],
        },
    }
    print(H.RESULT_MARK + json.dumps(result), flush=True)


# ===========================================================================
# 父进程编排
# ===========================================================================


class CoreWatcher(threading.Thread):
    """轮询发现本进程树下的 swanlab-core 子进程（累积合并，不覆写）。"""

    def __init__(self, parent_pid: int, expected: int, timeout: float) -> None:
        super().__init__(daemon=True, name="bench-core-watch")
        self.parent_pid = parent_pid
        self.expected = expected
        self.deadline = time.monotonic() + timeout
        self._lock = threading.Lock()
        self._pids: set[int] = set()
        self._halt = threading.Event()

    @property
    def pids(self) -> List[int]:
        with self._lock:
            return sorted(self._pids)

    def add(self, pids: List[int]) -> None:
        if not pids:
            return
        with self._lock:
            self._pids.update(int(p) for p in pids if p)

    def run(self) -> None:
        while not self._halt.is_set() and time.monotonic() < self.deadline:
            found = H.core_pids(parent_pid=self.parent_pid)
            if found:
                self.add(found)
            with self._lock:
                done = len(self._pids) >= self.expected
            if done:
                return
            self._halt.wait(0.05)

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=2)


def _child_argv(script: Path, args: argparse.Namespace, interpreter: str, core: str) -> List[str]:
    return [
        interpreter,
        str(script),
        "--worker-none" if core == "none" else "--worker",
        "--steps",
        str(args.steps),
        "--keys",
        str(args.keys),
        "--record-batch",
        str(args.record_batch),
        "--record-interval",
        str(args.record_interval),
        "--project",
        args.project,
        "--experiment",
        args.experiment,
        *(["--skip-store"] if args.skip_store else []),
        "--media-every",
        str(args.media_every),
        "--save-count",
        str(args.save_count),
    ]


def run_case(args: argparse.Namespace, core: str, workdir: Path, interpreter: str) -> Dict[str, Any]:
    """跑一次 case：起 mock → 起子进程 → 采样 → 校验 → 汇总。"""
    H.cleanup_stale_cores()
    case_dir = workdir / f"{core}-{args.scenario}-{int(time.time() * 1000) % 100000}"
    child_cwd = case_dir / "cwd"
    child_root = case_dir / "root"
    for path in (child_cwd, child_root, case_dir / "mock"):
        path.mkdir(parents=True, exist_ok=True)

    mock = H.MockCloud(case_dir / "mock", latency_ms=args.latency_ms)
    script = Path(__file__).resolve()
    checks = H.Checks()
    phases: Dict[str, float] = {}
    spawn_at = 0.0
    exit_at = 0.0
    sampler: Optional[H.Sampler] = None
    watcher: Optional[CoreWatcher] = None
    discovered_core_pids: set[int] = set()
    core_lock = threading.Lock()
    proc_pid: Optional[int] = None

    with mock:
        env = H.build_child_env(
            mock.url,
            root=child_root,
            cwd=child_cwd,
            core_impl=None if core == "none" else core,
        )
        stderr_path = case_dir / "child.stderr.log"
        stdout_lines: List[str] = []

        def pids_fn() -> List[int]:
            with core_lock:
                cores = set(discovered_core_pids)
            if watcher is not None:
                cores.update(watcher.pids)
            pids = [p for p in (proc_pid, *sorted(cores), mock.pid) if p]
            return pids

        sampler = H.Sampler(pids_fn, interval=0.05)
        sampler.start()
        spawn_at = sampler.capture()

        with stderr_path.open("w", encoding="utf-8") as stderr_file:
            proc = subprocess.Popen(
                _child_argv(script, args, interpreter, core),
                cwd=str(child_cwd),
                env=env,
                stdout=subprocess.PIPE,
                stderr=stderr_file,
                text=True,
                bufsize=1,
            )
            proc_pid = proc.pid
            expected_cores = 0 if core != "go" else 1
            if expected_cores:
                watcher = CoreWatcher(proc.pid, expected_cores, timeout=60.0)
                watcher.start()

            assert proc.stdout is not None
            for line in proc.stdout:
                stdout_lines.append(line)
                if line.startswith(H.PHASE_MARK):
                    payload = json.loads(line[len(H.PHASE_MARK) :])
                    reported_cores = [int(p) for p in (payload.get("core_pids") or []) if p]
                    if reported_cores:
                        with core_lock:
                            discovered_core_pids.update(reported_cores)
                        if watcher is not None:
                            watcher.add(reported_cores)
                    phases[str(payload.get("name"))] = sampler.capture()
            proc.wait(timeout=args.timeout)
            exit_at = time.monotonic()

        if watcher is not None:
            watcher.stop()
            with core_lock:
                discovered_core_pids.update(watcher.pids)
        if sampler is not None:
            sampler.stop()

        with core_lock:
            core_pids = sorted(discovered_core_pids)
        stderr_text = "".join(H.stderr_for_pid(pid) for pid in core_pids)
        child_stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
        mock_stats = mock.stats()

    result: Dict[str, Any] = {}
    for line in stdout_lines:
        if line.startswith(H.RESULT_MARK):
            result = json.loads(line[len(H.RESULT_MARK) :])

    if not result:
        raise RuntimeError(
            f"子进程未产出结果（core={core}, rc={proc.returncode}）\n"
            f"--- stdout ---\n{''.join(stdout_lines[-40:])}\n--- child stderr ---\n{child_stderr[-4000:]}\n"
            f"--- core stderr ---\n{stderr_text[-4000:]}"
        )

    # 在 sampler 停止之后由父进程离线扫描本地文件，避免回读操作污染被测进程 CPU/内存采样。
    if core != "none" and not args.skip_store and result.get("run_dir"):
        result["local_records"] = H.scan_local_records(Path(result["run_dir"]))

    # ---- 校验 -------------------------------------------------------------
    expect_impl = "none" if core == "none" else core
    checks.add("子进程未连真实后端", not result["isolation_problems"], str(result["isolation_problems"]))
    house_total = mock_stats["house_metrics_total"]
    extra_keys: List[str] = []
    if core == "go":
        checks.equal("Go core 进程数", len(core_pids), 1)
    elif core == "python":
        checks.equal("Python core 不应启动 core 进程", len(core_pids), 0)
    if core != "none":
        checks.equal("实际使用的 core 与预期一致", result["core_impl"], expect_impl)
        checks.equal("mock 收到标量条数", mock_stats["metrics_by_type"]["scalar"], result["scalar_expected"])
        local = result["local_records"]
        if not args.skip_store:
            checks.equal("本地文件标量条数", local.get("scalar", 0), result["scalar_expected"])
            # 这四类记录在两个 core 里都是一条落盘对应一条请求，因此本地条数应等于 mock
            # 收到的条数。column 走 /series、save 走文件或 profile，不在此列。
            for kind in ("scalar", "log", "media", "histogram"):
                checks.equal(
                    f"mock {kind} 条数与本地落盘一致",
                    mock_stats["metrics_by_type"].get(kind, 0),
                    local.get(kind, 0),
                )
            # uploaded_records 统计已处理完的记录：上报的四类 + column + save。
            # save 只要处理完就计数，和 mock 的 files_completed 不是一回事。
            delivered = house_total + local.get("column", 0) + local.get("save", 0)
            uploaded = result["upload"]["max_uploaded_records"]
            # uploaded_records 是每 0.2s 采样一次的快照，收尾阶段最后一刷可能晚于末次采样，
            # 因此只保证「不会超过实际交付的条数」；超出说明有记录处理完却没上传。
            checks.add("上传进度未超过已交付条数", uploaded <= delivered, f"uploaded={uploaded} delivered={delivered}")
            if args.media_every > 0:
                expected_media = len(range(0, int(args.steps), int(args.media_every)))
                checks.equal(
                    "mock media 条数与预期一致",
                    mock_stats["metrics_by_type"].get("media", 0),
                    expected_media,
                )
            if args.save_count > 0:
                # 用户 save 走 /files/prepare + /files/complete，每个文件一次 complete；
                # 配置/requirements/metadata 这类走 profile 接口，不调用 complete。
                checks.equal("用户 save 文件全部上传完成", mock_stats["files_completed"], int(args.save_count))
        # 列定义校验：mock 在 /series 收到的 key 必须覆盖 run.log 用过的每一个 key。
        expected_keys = {f"metric_{index}" for index in range(int(args.keys))}
        series_keys = set(mock_stats["series_keys"])
        missing_keys = sorted(expected_keys - series_keys)
        checks.add("series 覆盖全部 log key", not missing_keys, f"缺失 {missing_keys[:10]}")
        extra_keys = sorted(series_keys - expected_keys)
        if args.latency_ms == 0:
            checks.add("上传全程未被暂停", not result["upload"]["paused_seen"], "出现过 upload.paused=true")
        checks.add(
            "core 日志无报错",
            not any(token in stderr_text for token in ("panic", "store errored", "Error")),
            stderr_text[-800:],
        )
        checks.add("子进程无异常堆栈", "Traceback" not in child_stderr, child_stderr[-800:])

    # ---- 采样归并 ---------------------------------------------------------
    child_pid = proc.pid
    all_pids = [child_pid, *core_pids]
    mock_pids = [mock.pid] if mock.pid else []
    init_at = phases.get("init_done", spawn_at)
    finish_at = phases.get("finish_done", exit_at)
    wall_s = max(exit_at - spawn_at, 1e-9)
    cpu_window_s = max(finish_at - spawn_at, 1e-9)
    cpu_core_s = sampler.cpu_seconds(core_pids, spawn_at, finish_at) if sampler is not None else 0.0
    cpu_s = sampler.cpu_seconds(all_pids, spawn_at, finish_at) if sampler is not None else 0.0
    mock_cpu_s = sampler.cpu_seconds(mock_pids, spawn_at, finish_at) if sampler is not None else 0.0

    record = {
        "core": core,
        "scenario": args.scenario,
        "profile": args.profile,
        "skip_store": bool(args.skip_store),
        "steps": result["steps"],
        "keys": result["keys"],
        "record_batch": result["record_batch"],
        "record_interval": result["record_interval"],
        "scalar_expected": result["scalar_expected"],
        "media_logged": result.get("media_logged", 0),
        "save_files": result.get("save_files", 0),
        "init_ms": result["init_ms"],
        "producer_s": result["producer_s"],
        "finish_ms": result["finish_ms"],
        "total_s": result["total_s"],
        "wall_s": round(wall_s, 4),
        "producer_rec_per_s": round(result["scalar_expected"] / max(result["producer_s"], 1e-9), 1),
        "e2e_rec_per_s": round(result["scalar_expected"] / max(result["total_s"], 1e-9), 1),
        "latency_us": result["latency_us"],
        "latency_head_us": result["latency_head_us"],
        "has_phys_footprint": sampler.has_phys_footprint if sampler else False,
        "rss_child_peak": sampler.peak_rss([child_pid], init_at, finish_at) if sampler else 0,
        "rss_core_peak": sampler.peak_rss(core_pids, init_at, finish_at) if sampler else 0,
        "rss_total_peak": sampler.peak_rss(all_pids, spawn_at, finish_at) if sampler else 0,
        "phys_child_peak": sampler.peak_phys([child_pid], init_at, finish_at) if sampler else 0,
        "phys_core_peak": sampler.peak_phys(core_pids, init_at, finish_at) if sampler else 0,
        "phys_total_peak": sampler.peak_phys(all_pids, spawn_at, finish_at) if sampler else 0,
        "cpu_core_s": round(cpu_core_s, 2),
        "cpu_total_s": round(cpu_s, 2),
        "cpu_pct_total": round(100.0 * cpu_s / cpu_window_s, 1),
        "run_file_bytes": result["run_file_bytes"],
        "local_records": result["local_records"],
        "mock": {
            "by_type": mock_stats["metrics_by_type"],
            "house_total": house_total,
            "batches": mock_stats["metrics_batches"],
            "bytes": mock_stats["metrics_bytes"],
            "files_completed": mock_stats["files_completed"],
            "series_keys": len(mock_stats["series_keys"]),
            "series_requests": mock_stats["series_requests"],
            "series_extra_keys": extra_keys,
            "experiments": len(mock_stats["experiments"]),
            "heartbeats": mock_stats["heartbeats"],
            "stop_experiments": mock_stats["stop_experiments"],
            "cpu_s": round(mock_cpu_s, 2),
            "cpu_pct": round(100.0 * mock_cpu_s / cpu_window_s, 1),
            "rss_peak": sampler.peak_rss(mock_pids, spawn_at, finish_at) if sampler else 0,
            "phys_peak": sampler.peak_phys(mock_pids, spawn_at, finish_at) if sampler else 0,
        },
        "upload": result["upload"],
        "checks": [{"name": c.name, "ok": c.ok, "detail": c.detail} for c in checks.items],
        "ok": checks.ok(),
    }
    if not args.keep:
        shutil.rmtree(case_dir, ignore_errors=True)
    return record


# ===========================================================================
# 报告
# ===========================================================================


def _fmt_bytes(value: int) -> str:
    if value >= 1 << 30:
        return f"{value / (1 << 30):.2f} GB"
    if value >= 1 << 20:
        return f"{value / (1 << 20):.2f} MB"
    if value >= 1 << 10:
        return f"{value / (1 << 10):.1f} KB"
    return f"{value} B"


def render(records: List[Dict[str, Any]]) -> str:
    rows: List[List[Any]] = []
    for rec in records:
        label = rec["core"] + ("/skip_store" if rec["skip_store"] else "")
        rows.append(
            [
                label,
                rec["init_ms"],
                rec["producer_s"],
                rec["finish_ms"],
                rec["total_s"],
                rec["producer_rec_per_s"],
                rec["e2e_rec_per_s"],
                rec["latency_us"]["p99"],
                _fmt_bytes(rec["rss_total_peak"]),
                _fmt_bytes(rec.get("phys_total_peak", 0)) if rec.get("has_phys_footprint") else "—",
                _fmt_bytes(rec["run_file_bytes"]),
                f"{sum(1 for c in rec['checks'] if c['ok'])}/{len(rec['checks'])}",
            ]
        )
    headers = [
        "core",
        "init(ms)",
        "produce(s)",
        "finish(ms)",
        "total(s)",
        "prod rec/s",
        "e2e rec/s",
        "log p99(us)",
        "peak RSS",
        "phys RSS",
        "runfile",
        "checks",
    ]
    out = [H.render_table(headers, rows, widths=[16, 12, 12, 12, 11, 13, 13, 13, 12, 12, 11, 9])]
    out.append("")
    out.append("  上传状态（0.2s 采样）与 mock 计数：")
    detail_rows = []
    for rec in records:
        label = rec["core"] + ("/skip_store" if rec["skip_store"] else "")
        detail_rows.append(
            [
                label,
                rec["upload"]["max_buffered_bytes"],
                rec["upload"]["max_pending_records"],
                rec["upload"]["max_uploaded_records"],
                rec["mock"]["house_total"],
                rec["mock"]["by_type"].get("media", 0),
                rec["mock"]["files_completed"],
                rec["mock"]["batches"],
                "YES" if rec["upload"]["paused_seen"] else "no",
                rec["mock"]["series_requests"],
            ]
        )
    headers2 = [
        "core",
        "max buf(B)",
        "max pending",
        "uploaded(采样)",
        "mock house",
        "media",
        "files",
        "batches",
        "paused",
        "series",
    ]
    out.append(H.render_table(headers2, detail_rows, widths=[16, 12, 13, 12, 12, 7, 7, 9, 9, 8]))
    return "\n".join(out)


# ===========================================================================
# CLI
# ===========================================================================


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SwanLab L2 end-to-end benchmark (online + local mock cloud)")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker-none", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--record-batch", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--record-interval", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--scenario", choices=sorted(SCENARIOS), default="C")
    parser.add_argument("--steps", type=int, default=None, help="覆盖场景步数")
    parser.add_argument("--keys", type=int, default=None, help="覆盖场景指标数")
    parser.add_argument("--core", choices=("none", "python", "go", "both"), default="both")
    parser.add_argument("--profile", choices=sorted(PROFILES), default="default")
    parser.add_argument("--skip-store", action="store_true", help="对比 core.skip_store=True")
    parser.add_argument("--media-every", type=int, default=0, help="每 N 步记一张 16x16 图片（0=不开）")
    parser.add_argument("--save-count", type=int, default=0, help="run.save 的文件数（policy=end，0=不开）")
    parser.add_argument("--repeats", type=int, default=1, help="每个 case 重复次数，取中位数")
    parser.add_argument("--latency-ms", type=int, default=0, help="mock 响应注入的固定延迟")
    parser.add_argument("--interpreter", default=None, help="被测子进程解释器（用于抽跑旧版本）")
    parser.add_argument("--project", default="bench-project")
    parser.add_argument("--experiment", default="bench")
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--keep", action="store_true", help="保留子进程工作目录")
    parser.add_argument("--json", default=None, help="结果落盘文件名（写入 tests/benchmark/results）")
    args = parser.parse_args(argv)

    preset = SCENARIOS[args.scenario]
    args.steps = args.steps or preset["steps"]
    args.keys = args.keys or preset["keys"]
    profile = PROFILES[args.profile]
    if args.record_batch is None:
        args.record_batch = int(profile["record_batch"])
    if args.record_interval is None:
        args.record_interval = float(profile["record_interval"])
    return args


def _median(values: List[float]) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    return ordered[len(ordered) // 2]


def pick_best(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """同一 case 多次重复：取 total_s 中位数那一次作为代表。"""
    if len(records) == 1:
        return records[0]
    target = _median([rec["total_s"] for rec in records])
    return min(records, key=lambda rec: abs(rec["total_s"] - target))


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.worker_none:
        run_worker_none(args)
        return 0
    if args.worker:
        run_worker(args)
        return 0

    interpreter = args.interpreter or sys.executable
    cores = CORES if args.core == "both" else (args.core,)
    workdir = Path(tempfile.mkdtemp(prefix="swanlab-bench-e2e-"))
    print(f"  SwanLab L2 e2e benchmark  scenario={args.scenario} steps={args.steps} keys={args.keys}")
    print(
        f"  profile={args.profile} record_batch={args.record_batch} "
        f"record_interval={args.record_interval} skip_store={args.skip_store} repeats={args.repeats}"
    )
    print(f"  core={','.join(cores)}  workdir={workdir}")

    environment = H.environment_info()
    records: List[Dict[str, Any]] = []
    try:
        for core in cores:
            runs: List[Dict[str, Any]] = []
            for index in range(max(1, args.repeats)):
                print(f"  ... core={core} repeat={index + 1}/{args.repeats}", flush=True)
                runs.append(run_case(args, core, workdir, interpreter))
            best = pick_best(runs)
            best["all_totals_s"] = [rec["total_s"] for rec in runs]
            records.append(best)
    finally:
        if not args.keep:
            shutil.rmtree(workdir, ignore_errors=True)

    print()
    print(render(records))
    print()
    for rec in records:
        failed = [c for c in rec["checks"] if not c["ok"]]
        label = rec["core"] + ("/skip_store" if rec["skip_store"] else "")
        if failed:
            print(f"  ✗ {label}: {len(failed)}/{len(rec['checks'])} 项校验未通过")
            for item in failed:
                print(f"    - {item['name']}: {item['detail']}")
        else:
            print(f"  ✓ {label}: 校验全部通过（{len(rec['checks'])} 项）")

    payload = {"env": environment, "args": vars(args), "records": records}
    name = args.json or f"e2e-{args.scenario}-{args.profile}"
    path = H.write_result(payload, name)
    print(f"\n  结果已写入 {path}")

    return 0 if all(rec["ok"] for rec in records) else 1


if __name__ == "__main__":
    sys.exit(main())
