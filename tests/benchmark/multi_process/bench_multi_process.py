"""
@author: nexisato
@file: bench_multi_process.py
@description: 多进程并行基准：N 个 SDK 进程写同一个实验（parallel="shared"），对比
              Go 共享模式 / Go 独占模式 / Python core 三种运行方式。

三组用来分离两件事：
  - core 实现差异：Go 独占（N 个 core 进程） vs Python core（每进程各自实现）
  - 运行方式差异：Go 共享（1 个 core 进程） vs Go 独占（N 个 core 进程）

注意 detach=True 只共享 Go 进程、HTTP 连接池和站点探测结果；本地写入、上传队列、
预算、心跳、步数去重都是每个会话各一份，所以 N 个 rank 是「1 个进程 + N 套会话资源」。

用法::

    uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 4 --steps 20 --keys 200
    uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 16 --steps 100 --keys 1000
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import common as H  # noqa: E402

# 三种运行方式；go-shared 用 detach=True，go-owner 用默认的每进程独占。
TOPOLOGIES: Dict[str, Dict[str, Any]] = {
    "go-shared": {"core": "go", "detach": True, "expected_cores": "shared"},
    "go-owner": {"core": "go", "detach": False, "expected_cores": "per-rank"},
    "python": {"core": "python", "detach": False, "expected_cores": 0},
}

DEFAULT_STEPS = 100
DEFAULT_KEYS = 1000
DEFAULT_RANKS = 4


# ===========================================================================
# 子进程 Worker
# ===========================================================================


def _emit_phase(name: str, **extra: Any) -> None:
    payload = {"name": name, "t": time.monotonic(), "rank": extra.pop("rank", -1)}
    payload.update(extra)
    print(H.PHASE_MARK + json.dumps(payload), flush=True)


def _cores_seen(detach: bool) -> List[int]:
    """返回本进程实际连着的 core pid。

    共享模式下 core 不是本进程的子进程，直接读共享 port-file；
    独占模式下 core 是子进程，按 ppid 反查。
    """
    if detach:
        port = Path(f"/tmp/swanlab-core-{os.getuid()}") / "core.port"
        pid = H.read_port_file_pid(port)
        return [pid] if pid else []
    return H.core_pids(parent_pid=os.getpid())


def run_worker(args: argparse.Namespace) -> None:
    import swanlab

    core_impl = H.backend_kind()
    violations = H.isolation_problems()
    _emit_phase("ready", rank=args.rank, core_impl=core_impl)

    settings = swanlab.Settings(
        core=swanlab.Settings.Core(
            detach=bool(args.detach),
            record_batch=int(args.record_batch),
            record_interval=float(args.record_interval),
        ),
        probe=swanlab.Settings.Probe(monitor=False),
    )

    t_init_start = time.perf_counter()
    run = swanlab.init(
        mode="online",
        project=args.project,
        experiment_name=args.experiment,
        id=args.run_id,
        parallel="shared",
        settings=settings,
    )
    init_ms = (time.perf_counter() - t_init_start) * 1000.0
    _emit_phase("init_done", rank=args.rank, init_ms=round(init_ms, 2))
    cores_seen = _cores_seen(bool(args.detach))

    probe = H.StatsProbe(run)
    probe.start()

    combine_keys = [f"combine/train_{index}" for index in range(int(args.keys))]
    independent_keys = [f"independent_{args.rank}/train_{index}" for index in range(int(args.keys))]
    ranks = int(args.ranks)
    latencies: List[float] = []
    error: Optional[str] = None

    try:
        t_total_start = time.perf_counter()
        for step in range(int(args.steps)):
            combine = {key: 0.5 + (step % 1000) * 1e-4 for key in combine_keys}
            independent = {key: 0.25 + (step % 1000) * 1e-4 for key in independent_keys}
            t0 = time.perf_counter()
            run.log(combine, step=step * ranks + args.rank)
            run.log(independent, step=step)
            latencies.append((time.perf_counter() - t0) * 1e6)
            if args.sleep:
                time.sleep(float(args.sleep))
            if step % 20 == 0:
                _emit_phase("progress", rank=args.rank, step=step)

        producer_s = time.perf_counter() - t_total_start
        _emit_phase("producer_done", rank=args.rank, producer_s=round(producer_s, 4))

        t_finish_start = time.perf_counter()
        run.finish()
        finish_ms = (time.perf_counter() - t_finish_start) * 1000.0
        total_s = time.perf_counter() - t_total_start
        _emit_phase("finish_done", rank=args.rank, finish_ms=round(finish_ms, 2))

        probe.stop()
        run_dir = Path(run._ctx.run_dir)
        local_counts = H.scan_local_records(run_dir)
        run_file_bytes = sum(path.stat().st_size for path in run_dir.glob("run-*.swanlab"))
        run_file = next(iter(sorted(run_dir.glob("run-*.swanlab")))).name if run_file_bytes else ""
    except BaseException as exc:  # noqa: BLE001 - 带上下文回传父进程
        error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        probe.stop()

    result = {
        "ok": error is None,
        "error": error,
        "rank": int(args.rank),
        "ranks": ranks,
        "steps": int(args.steps),
        "keys": int(args.keys),
        "sleep": float(args.sleep),
        "detach": bool(args.detach),
        "core_impl": core_impl,
        "isolation_problems": violations,
        "init_ms": round(init_ms, 2),
        "producer_s": round(producer_s, 4),
        "finish_ms": round(finish_ms, 2),
        "total_s": round(total_s, 4),
        "scalar_expected": int(args.steps) * int(args.keys) * 2,
        "latency_us": H.latency_stats(latencies),
        "local_records": local_counts,
        "run_file": run_file,
        "run_file_bytes": run_file_bytes,
        "upload": probe.summary(),
        "cores_seen": cores_seen,
        "cores_seen_at_exit": _cores_seen(bool(args.detach)),
    }
    print(H.RESULT_MARK + json.dumps(result), flush=True)


# ===========================================================================
# 父进程编排
# ===========================================================================


def _child_argv(script: Path, args: argparse.Namespace, rank: int, interpreter: str, detach: bool) -> List[str]:
    return [
        interpreter,
        str(script),
        "--worker",
        "--rank",
        str(rank),
        "--ranks",
        str(args.ranks),
        "--steps",
        str(args.steps),
        "--keys",
        str(args.keys),
        "--sleep",
        str(args.sleep),
        "--run-id",
        args.run_id,
        "--project",
        args.project,
        "--experiment",
        args.experiment,
        "--record-batch",
        str(args.record_batch),
        "--record-interval",
        str(args.record_interval),
        *(["--detach"] if detach else []),
    ]


def run_case(args: argparse.Namespace, topology: str, workdir: Path, interpreter: str) -> Dict[str, Any]:
    spec = TOPOLOGIES[topology]
    H.cleanup_stale_cores()

    case_dir = workdir / f"{topology}-{int(time.time() * 1000) % 100000}"
    mock_dir = case_dir / "mock"
    mock_dir.mkdir(parents=True, exist_ok=True)
    checks = H.Checks()
    script = Path(__file__).resolve()
    procs: List[subprocess.Popen] = []
    stderr_paths: List[Path] = []
    stdout_lines: List[str] = []
    phases: Dict[str, float] = {}
    spawn_at = 0.0
    exit_at = 0.0
    core_pids: List[int] = []
    core_ever_seen: set = set()

    mock = H.MockCloud(mock_dir, latency_ms=args.latency_ms)
    with mock:
        spawn_at = time.monotonic()
        for rank in range(args.ranks):
            rank_dir = case_dir / f"rank{rank}"
            (rank_dir / "cwd").mkdir(parents=True, exist_ok=True)
            (rank_dir / "root").mkdir(parents=True, exist_ok=True)
            env = H.build_child_env(
                mock.url,
                root=rank_dir / "root",
                cwd=rank_dir / "cwd",
                core_impl=spec["core"],
            )
            stderr_path = rank_dir / "child.stderr.log"
            stderr_paths.append(stderr_path)
            handle = stderr_path.open("w", encoding="utf-8")
            procs.append(
                subprocess.Popen(
                    _child_argv(script, args, rank, interpreter, bool(spec["detach"])),
                    cwd=str(rank_dir / "cwd"),
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=handle,
                    text=True,
                    bufsize=1,
                )
            )
            # 让 rank 依次启动，避免同时抢启动锁（真实 DDP 也是逐个拉起）。
            time.sleep(0.15)

        sampler: Optional[H.Sampler] = H.Sampler(lambda: [p.pid for p in procs] + core_pids)
        sampler.start()

        # 汇总所有 rank 的 stdout（并发读取，避免管道互锁）
        import threading

        def _drain(proc: subprocess.Popen) -> None:
            assert proc.stdout is not None
            for line in proc.stdout:
                stdout_lines.append(line)
                if line.startswith(H.PHASE_MARK):
                    payload = json.loads(line[len(H.PHASE_MARK) :])
                    phases[f"{payload.get('rank')}:{payload.get('name')}"] = time.monotonic()

        readers = [threading.Thread(target=_drain, args=(proc,), daemon=True) for proc in procs]
        for reader in readers:
            reader.start()

        # 轮询发现 core 进程
        deadline = time.monotonic() + 120.0
        while time.monotonic() < deadline:
            found = H.core_pids()
            if found:
                core_pids = found
                core_ever_seen.update(found)
            if all(proc.poll() is not None for proc in procs):
                break
            time.sleep(0.1)

        for proc in procs:
            proc.wait(timeout=args.timeout)
        exit_at = time.monotonic()
        for reader in readers:
            reader.join(timeout=10)
        sampler.stop()
        final_core_pids = H.core_pids()
        core_ever_seen.update(final_core_pids)
        stderr_text = "".join(H.stderr_for_pid(pid) for pid in sorted(core_ever_seen))
        mock_stats = mock.stats()

    for proc in procs:
        if proc.stderr is not None:
            proc.stderr.close()

    results: List[Dict[str, Any]] = []
    for line in stdout_lines:
        if line.startswith(H.RESULT_MARK):
            results.append(json.loads(line[len(H.RESULT_MARK) :]))
    results.sort(key=lambda item: item["rank"])

    if len(results) != args.ranks:
        child_stderr = "\n".join(path.read_text(encoding="utf-8", errors="replace")[-1500:] for path in stderr_paths)
        raise RuntimeError(
            f"仅 {len(results)}/{args.ranks} 个 rank 产出结果（topology={topology}）\n"
            f"--- child stderr ---\n{child_stderr}\n--- core stderr ---\n{stderr_text[-3000:]}"
        )

    # ---- 校验 -------------------------------------------------------------
    expected_scalars = args.ranks * args.steps * args.keys * 2
    local_scalars = sum(item["local_records"].get("scalar", 0) for item in results)
    mock_scalars = mock_stats["metrics_by_type"]["scalar"]

    if spec["core"] != "none":
        checks.add(
            "实际使用的 core 与预期一致",
            all(item["core_impl"] == spec["core"] for item in results),
            str(sorted({item["core_impl"] for item in results})),
        )
        checks.add(
            "子进程未连真实后端",
            all(not item["isolation_problems"] for item in results),
            str([item["isolation_problems"] for item in results if item["isolation_problems"]]),
        )
    # 判断依据是子进程自己上报的 pid；父进程轮询可能错过短命进程。
    reported: List[int] = []
    for item in results:
        reported.extend(item.get("cores_seen") or [])
    if spec["expected_cores"] == "shared":
        checks.equal("共享模式下只有 1 个 core 进程", len(set(reported)), 1)
        checks.add("所有 rank 连的是同一个 core", len(set(reported)) == 1, f"pids={sorted(set(reported))}")
    elif spec["expected_cores"] == "per-rank":
        checks.equal("每个 rank 各启动 1 个 core 进程", len(reported), args.ranks)
        checks.equal("各 rank 用的 core 进程互不相同", len(set(reported)), args.ranks)
    else:
        checks.equal("Python core 不应启动 core 进程", len(set(reported)), 0)
    checks.equal("mock 收到标量条数", mock_scalars, expected_scalars)
    checks.equal("本地落盘标量条数合计", local_scalars, expected_scalars)
    checks.add(
        "各 rank 均无上传暂停",
        not any(item["upload"]["paused_seen"] for item in results),
        "出现 upload.paused=true",
    )
    checks.add(
        "core stderr 无错误",
        not any(token in stderr_text for token in ("panic", "store errored", "Error")),
        stderr_text[-800:],
    )
    checks.add(
        "各 rank 落盘条数与预期一致",
        _local_counts_ok(results, args),
        f"每 rank 期望 {args.steps * args.keys * 2} 条标量",
    )

    all_pids = [proc.pid for proc in procs] + core_pids
    child_pids = [proc.pid for proc in procs]
    wall_s = max(exit_at - spawn_at, 1e-9)
    cpu_s = sampler.cpu_seconds(all_pids)
    p99 = max(item["latency_us"]["p99"] for item in results)

    record = {
        "topology": topology,
        "core": spec["core"],
        "detach": bool(spec["detach"]),
        "ranks": args.ranks,
        "steps": args.steps,
        "keys": args.keys,
        "sleep": float(args.sleep),
        "record_batch": int(args.record_batch),
        "record_interval": float(args.record_interval),
        "wall_s": round(wall_s, 4),
        "scalar_expected": expected_scalars,
        "throughput_rec_per_s": round(expected_scalars / wall_s, 1),
        "init_ms_median": _median([item["init_ms"] for item in results]),
        "producer_s_max": max(item["producer_s"] for item in results),
        "finish_ms_max": max(item["finish_ms"] for item in results),
        "latency_p99_max_us": p99,
        "rss_children_peak": sampler.peak_rss(child_pids),
        "rss_core_peak": sampler.peak_rss(core_pids),
        "rss_total_peak": sampler.peak_rss(all_pids),
        "cpu_core_s": round(sampler.cpu_seconds(core_pids), 2),
        "cpu_total_s": round(cpu_s, 2),
        "cpu_pct_total": round(100.0 * cpu_s / wall_s, 1),
        "core_count": len(set(reported)),
        "core_pids": sorted(set(reported)),
        "disk_bytes": sum(item["run_file_bytes"] for item in results),
        "mock": {
            "by_type": mock_stats["metrics_by_type"],
            "house_total": mock_stats["house_metrics_total"],
            "batches": mock_stats["metrics_batches"],
            "experiments": len(mock_stats["experiments"]),
            "heartbeats": mock_stats["heartbeats"],
            "files_completed": mock_stats["files_completed"],
        },
        "per_rank": [
            {
                "rank": item["rank"],
                "init_ms": item["init_ms"],
                "producer_s": item["producer_s"],
                "finish_ms": item["finish_ms"],
                "p99_us": item["latency_us"]["p99"],
                "scalar_local": item["local_records"].get("scalar", 0),
                "uploaded_records": item["upload"]["max_uploaded_records"],
                "max_buffered_bytes": item["upload"]["max_buffered_bytes"],
                "max_pending_records": item["upload"]["max_pending_records"],
                "paused": item["upload"]["paused_seen"],
            }
            for item in results
        ],
        "checks": [{"name": c.name, "ok": c.ok, "detail": c.detail} for c in checks.items],
        "ok": checks.ok(),
    }
    H.cleanup_stale_cores()
    if not args.keep:
        shutil.rmtree(case_dir, ignore_errors=True)
    return record


def _local_counts_ok(results: List[Dict[str, Any]], args: argparse.Namespace) -> bool:
    """每个 rank 每步写两类指标各 keys 条，本地标量条数应为 steps * keys * 2。

    共享指标用 step * ranks + rank 生成，各 rank 的 step 互不重叠，因此只要条数
    对得上就不会互相覆盖。
    """
    expected = args.steps * args.keys * 2
    return expected > 0 and all(item["local_records"].get("scalar", 0) == expected for item in results)


# ===========================================================================
# 报告
# ===========================================================================


def _median(values: List[float]) -> float:
    ordered = sorted(values)
    return ordered[len(ordered) // 2] if ordered else 0.0


def _spread(values: List[float]) -> float:
    """同一配置多次重复的相对离散度：(max-min)/中位数，百分比。"""
    if len(values) < 2:
        return 0.0
    middle = _median(values)
    return round((max(values) - min(values)) / middle * 100, 1) if middle else 0.0


def pick_best(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """同一拓扑多次重复：取墙钟中位数那一次作为代表。"""
    if len(records) == 1:
        return records[0]
    target = _median([rec["wall_s"] for rec in records])
    return min(records, key=lambda rec: abs(rec["wall_s"] - target))


def _fmt_bytes(value: int) -> str:
    for unit, scale in (("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)):
        if value >= scale:
            return f"{value / scale:.1f} {unit}"
    return f"{value} B"


def render(records: List[Dict[str, Any]]) -> str:
    rows = []
    for rec in records:
        rows.append(
            [
                rec["topology"],
                rec["ranks"],
                rec["core_count"],
                rec["wall_s"],
                f"{rec['spread_pct']}%",
                rec["throughput_rec_per_s"],
                rec["init_ms_median"],
                rec["finish_ms_max"],
                rec["latency_p99_max_us"],
                _fmt_bytes(rec["rss_children_peak"]),
                _fmt_bytes(rec["rss_core_peak"]),
                _fmt_bytes(rec["rss_total_peak"]),
                f"{rec['cpu_pct_total']}%",
                f"{sum(1 for c in rec['checks'] if c['ok'])}/{len(rec['checks'])}",
            ]
        )
    headers = [
        "topology",
        "ranks",
        "cores",
        "wall(s)",
        "spread",
        "rec/s",
        "init(ms)",
        "finish(ms)",
        "p99(us)",
        "rss py",
        "rss core",
        "rss total",
        "cpu",
        "checks",
    ]
    out = [
        H.render_table(
            headers,
            rows,
            widths=[14, 7, 7, 12, 8, 12, 10, 12, 11, 10, 10, 11, 8, 9],
        )
    ]
    out.append("")
    detail = []
    for rec in records:
        detail.append(
            [
                rec["topology"],
                rec["mock"]["by_type"]["scalar"],
                rec["mock"]["house_total"],
                rec["mock"]["batches"],
                rec["mock"]["experiments"],
                rec["mock"]["heartbeats"],
                _fmt_bytes(rec["disk_bytes"]),
            ]
        )
    out.append(
        H.render_table(
            ["topology", "mock scalar", "mock total", "batches", "experiments", "heartbeats", "disk"],
            detail,
            widths=[14, 13, 12, 10, 13, 12, 11],
        )
    )
    # per-rank 上传队列深度：max 为该拓扑各 rank 的峰值，用来判断 p99 抬升是否来自队列满。
    queue_rows = []
    for rec in records:
        ranks = rec.get("per_rank") or []
        queue_rows.append(
            [
                rec["topology"],
                rec["ranks"],
                max((item["max_buffered_bytes"] for item in ranks), default=0),
                max((item["max_pending_records"] for item in ranks), default=0),
                max((item["uploaded_records"] for item in ranks), default=0),
                sum(1 for item in ranks if item["paused"]),
            ]
        )
    out.append("")
    out.append(
        H.render_table(
            ["topology", "ranks", "max buf(B)", "max pending", "uploaded max", "paused ranks"],
            queue_rows,
            widths=[14, 7, 13, 13, 14, 13],
        )
    )
    return "\n".join(out)


# ===========================================================================
# CLI
# ===========================================================================


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SwanLab L3 multi-process benchmark (shared experiment)")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--rank", type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument("--detach", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    parser.add_argument("--keys", type=int, default=DEFAULT_KEYS, help="每步每类指标数（共享/独立各一份）")
    parser.add_argument("--sleep", type=float, default=0.0, help="每步 sleep 秒数；0.2 贴近真实训练")
    parser.add_argument("--topology", choices=sorted(TOPOLOGIES) + ["all"], default="all")
    parser.add_argument("--record-batch", type=int, default=10_000, help=argparse.SUPPRESS)
    parser.add_argument("--record-interval", type=float, default=1.5, help=argparse.SUPPRESS)
    parser.add_argument("--run-id", default=None, help="共享实验 id，默认自动生成")
    parser.add_argument("--project", default="bench-project")
    parser.add_argument("--experiment", default="bench-mp")
    parser.add_argument("--latency-ms", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=1, help="每个拓扑重复次数，取墙钟中位数")
    parser.add_argument("--interpreter", default=None)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--json", default=None)
    args = parser.parse_args(argv)
    args.run_id = args.run_id or f"bench-mp-{int(time.time())}"
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.worker:
        run_worker(args)
        return 0

    interpreter = args.interpreter or sys.executable
    topologies = sorted(TOPOLOGIES) if args.topology == "all" else (args.topology,)
    workdir = Path(tempfile.mkdtemp(prefix="swanlab-bench-mp-"))
    print(
        f"  SwanLab L3 multi-process benchmark  ranks={args.ranks} steps={args.steps} "
        f"keys={args.keys} sleep={args.sleep} run_id={args.run_id}"
    )
    print(f"  topologies={','.join(topologies)}  workdir={workdir}")

    environment = H.environment_info()
    records: List[Dict[str, Any]] = []
    try:
        for topology in topologies:
            runs: List[Dict[str, Any]] = []
            for index in range(max(1, args.repeats)):
                print(f"  ... topology={topology} repeat={index + 1}/{args.repeats}", flush=True)
                runs.append(run_case(args, topology, workdir, interpreter))
            best = pick_best(runs)
            best["all_wall_s"] = [rec["wall_s"] for rec in runs]
            best["spread_pct"] = _spread([rec["wall_s"] for rec in runs])
            records.append(best)
    finally:
        if not args.keep:
            shutil.rmtree(workdir, ignore_errors=True)

    print()
    print(render(records))
    print()
    for rec in records:
        failed = [c for c in rec["checks"] if not c["ok"]]
        if failed:
            print(f"  ✗ {rec['topology']}: {len(failed)}/{len(rec['checks'])} 项校验未通过")
            for item in failed:
                print(f"    - {item['name']}: {item['detail']}")
        else:
            print(f"  ✓ {rec['topology']}: 校验全部通过（{len(rec['checks'])} 项）")

    payload = {"env": environment, "args": vars(args), "records": records}
    name = args.json or f"mp-{args.ranks}rank-sleep{args.sleep}"
    path = H.write_result(payload, name)
    print(f"\n  结果已写入 {path}")
    return 0 if all(rec["ok"] for rec in records) else 1


if __name__ == "__main__":
    sys.exit(main())
