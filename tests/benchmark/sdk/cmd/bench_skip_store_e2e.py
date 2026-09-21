"""
@author: nexisato
@file: bench_skip_store_e2e.py
@time: 2026/9/11
@description: core.skip_store 端到端基准测试（mock HTTP，对比吞吐量与延时）

在 online 模式下完整跑一遍 ``swanlab.init → log/log_image/save → finish``，
所有 HTTP 端点用 responses mock，对比 ``core.skip_store`` 开关对**用户线程延时**
与**端到端吞吐量**的影响。

被测对象是真实运行时链路：事件总线 → BackgroundConsumer → CorePython → 本地 store
+ Transport → HTTP 上传，而非单独的存储模块。

工作负载（每个场景一致）：
  - 标量：STEPS 步 × KEYS 个 key，逐 step 调用 ``run.log``
  - 媒体：每隔 IMAGE_EVERY 步调用一次 ``run.log_image``（小图），共约 STEPS/IMAGE_EVERY 张
  - 文件保存：SAVE_COUNT 个 CUSTOM save，policy=end（生成阶段只登记，finish 时上传）

关注指标：
  1. 生产阶段墙钟耗时与吞吐量（records/s）
  2. finish 排空 + 上传耗时 (ms)
  3. 端到端墙钟耗时 (s) 与吞吐量 (records/s)
  4. run.log 主线程单次调用延迟 mean / p50 / p95 / p99 (μs)
  5. 本地产物（run-*.swanlab 字节数、media 文件数、files 软链接数）

说明：run.log_image / run.save 的单次延时不做对比——它们会被后台
BackgroundConsumer / Transport 的 GIL 竞争严重污染，无法反映用户线程成本；
媒体与 save 的收益通过 finish 排空耗时与总墙钟体现。为使 producer / finish
两阶段边界确定，``record_interval`` 设为很大值，上传统一发生在 finish。

每个场景在独立子进程中运行（responses mock + 全局单例无法安全复用），
父进程取多次运行的最优值后汇总对比。

用法：
  uv run pytest tests/benchmark/sdk/cmd/bench_skip_store_e2e.py -v -s

  # 独立脚本（可调参数）
  uv run python tests/benchmark/sdk/cmd/bench_skip_store_e2e.py --steps 500 --keys 20
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional

RESULT_MARK = "##BENCH_JSON## "

# ===========================================================================
# Mock HTTP 常量（与 tests/unit/sdk/cmd/init/test_init_e2e.py 保持一致）
# ===========================================================================

API_HOST = "https://api.fake.swanlab.cn"
WEB_HOST = "https://test.swanlab.cn"
USERNAME = "test-user"
PROJECT = "test-project"
RUN_ID = "test-run-id"
EXPERIMENT_CUID = "test-experiment-cuid"
API_KEY = "test-api-key"

# ===========================================================================
# 默认负载参数
# ===========================================================================

STEPS = 1000
KEYS = 100
IMAGE_EVERY = 4
IMAGE_SIZE = 64  # 64x64x3 uint8
SAVE_COUNT = 50
SAVE_SIZE = 4 * 1024

# 让 Transport 不在生产阶段中途 drain，保证 producer / finish 两阶段边界确定，
# 生产阶段只体现“事件入队 + 消费线程本地 store”的成本，finish 阶段统一上传。
RECORD_INTERVAL = 3600.0
REPEATS = 2


# ===========================================================================
# 统计辅助
# ===========================================================================


def percentile(sorted_data: List[float], pct: float) -> float:
    if not sorted_data:
        return 0.0
    idx = int(len(sorted_data) * pct / 100.0)
    return sorted_data[min(idx, len(sorted_data) - 1)]


def _latency_stats(latencies_us: List[float]) -> Dict[str, float]:
    if not latencies_us:
        return {"mean": 0.0, "p50": 0.0, "p95": 0.0, "p99": 0.0}
    s = sorted(latencies_us)
    return {
        "mean": round(sum(s) / len(s), 2),
        "p50": round(percentile(s, 50), 2),
        "p95": round(percentile(s, 95), 2),
        "p99": round(percentile(s, 99), 2),
    }


# ===========================================================================
# HTTP Mock（动态预签名）
# ===========================================================================


def _install_http_mocks(rsps) -> None:
    """注册 init + 指标 + 媒体 + save 全流程端点，预签名 URL 按请求动态返回。"""
    import responses as responses_lib

    def ok(**extra):
        body = {"message": "ok"}
        body.update(extra)
        return (200, {"Content-Type": "application/json"}, json.dumps(body))

    def presigned_callback(request):
        payload = json.loads(request.body)
        urls = [f"https://storage.fake.swanlab.cn/media/{i}" for i in range(len(payload["paths"]))]
        return (200, {"Content-Type": "application/json"}, json.dumps({"urls": urls}))

    def prepare_callback(request):
        payload = json.loads(request.body)
        count = len(payload.get("files", []))
        urls = [f"https://storage.fake.swanlab.cn/save/{i}" for i in range(count)]
        return (200, {"Content-Type": "application/json"}, json.dumps({"urls": urls}))

    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/login/api_key",
        json={"sid": "mock-sid", "expiredAt": "2099-12-31T23:59:59.000Z", "userInfo": {"username": USERNAME}},
        status=200,
    )
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/projects/{USERNAME}",
        json={"name": PROJECT, "username": USERNAME, "path": f"/{USERNAME}/{PROJECT}"},
        status=201,
    )
    rsps.add(
        responses_lib.GET,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}",
        json={
            "cuid": "test-project-cuid",
            "name": PROJECT,
            "version": 1,
            "group": {"username": USERNAME},
            "username": USERNAME,
            "path": f"/{USERNAME}/{PROJECT}",
            "visibility": "PRIVATE",
            "_count": {"experiments": 0, "contributors": 1, "collaborators": 0, "clones": 0},
        },
        status=200,
    )
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/experiment",
        json={"cuid": EXPERIMENT_CUID, "slug": RUN_ID, "name": "test-experiment"},
        status=201,
    )
    rsps.add(
        responses_lib.PUT,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/runs/{EXPERIMENT_CUID}/state",
        json={"message": "ok"},
        status=200,
    )
    rsps.add(
        responses_lib.PUT,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/runs/{EXPERIMENT_CUID}/profile",
        json={"message": "ok"},
        status=200,
    )
    rsps.add(
        responses_lib.POST, f"{API_HOST}/api/projects/{USERNAME}/{PROJECT}/series", json={"message": "ok"}, status=200
    )
    rsps.add(responses_lib.POST, f"{API_HOST}/api/house/metrics", json={"message": "ok"}, status=200)
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/house/experiments/{EXPERIMENT_CUID}/heartbeat",
        json={"message": "ok"},
        status=200,
    )
    # 媒体：动态预签名 + 对象存储 PUT
    rsps.add_callback(responses_lib.POST, f"{API_HOST}/api/resources/presigned/put", callback=presigned_callback)
    rsps.add(responses_lib.PUT, re.compile(r"https://storage\.fake\.swanlab\.cn/media/.*"), body="", status=200)
    # 文件保存：动态 prepare + 对象存储 PUT + complete
    rsps.add_callback(
        responses_lib.POST, f"{API_HOST}/api/experiment/{EXPERIMENT_CUID}/files/prepare", callback=prepare_callback
    )
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/experiment/{EXPERIMENT_CUID}/files/complete",
        json={"message": "ok"},
        status=201,
    )
    rsps.add(responses_lib.PUT, re.compile(r"https://storage\.fake\.swanlab\.cn/save/.*"), body="", status=200)


# ===========================================================================
# 子进程 Worker
# ===========================================================================


def run_worker(skip_store: bool, steps: int, keys: int, save_count: int) -> None:
    import numpy as np
    import responses as responses_lib

    import swanlab
    from swanlab.sdk.cmd.login import login_raw
    from swanlab.sdk.cmd.merge_settings import merge_settings

    tmp_dir = tempfile.mkdtemp(prefix="swanlab_e2e_bench_")
    src_dir = os.path.join(tmp_dir, "src")
    os.makedirs(src_dir, exist_ok=True)
    save_files = []
    for i in range(save_count):
        path = os.path.join(src_dir, f"ckpt_{i:03d}.bin")
        with open(path, "wb") as f:
            f.write(bytes([i % 256]) * SAVE_SIZE)
        save_files.append(path)

    finish_ms = 0.0
    total_s = 0.0
    run_file_bytes = 0
    media_files = 0
    save_links = 0
    result: Dict[str, Any] = {}

    try:
        with responses_lib.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            _install_http_mocks(rsps)
            merge_settings({"api_host": API_HOST, "web_host": WEB_HOST, "probe": {"monitor": False}})
            login_raw(api_key=API_KEY, host=API_HOST, save=False, print_welcome=False)

            settings = swanlab.Settings(
                core=swanlab.Settings.Core(skip_store=skip_store, record_interval=RECORD_INTERVAL)
            )
            run = swanlab.init(mode="online", project=PROJECT, log_dir=tmp_dir, settings=settings)

            log_lat: List[float] = []
            num_image = 0

            t_total_start = time.perf_counter()

            # ---- 标量 + 媒体生产阶段 ----
            for s in range(steps):
                data = {f"k{i}": 0.5 + (s % 1000) * 1e-4 for i in range(keys)}
                t0 = time.perf_counter()
                run.log(data, step=s)
                log_lat.append((time.perf_counter() - t0) * 1e6)

                if s % IMAGE_EVERY == 0:
                    img = np.zeros((IMAGE_SIZE, IMAGE_SIZE, 3), dtype=np.uint8)
                    run.log_image(key="img", data=img, step=s)
                    num_image += 1

            t_producer_end = time.perf_counter()

            # ---- 文件保存（policy=end，仅登记，finish 上传）----
            for path in save_files:
                run.save(path, base_path=src_dir, policy="end")

            # ---- finish：排空 consumer + transport ----
            t_finish_start = time.perf_counter()
            run.finish()
            finish_ms = (time.perf_counter() - t_finish_start) * 1e6 / 1e3

            t_total_end = time.perf_counter()

            # 本地产物统计（非 skip 才有）
            run_dir = run._ctx.run_dir
            for p in run_dir.glob("run-*.swanlab"):
                run_file_bytes += p.stat().st_size
            media_dir = run_dir / "media"
            if media_dir.exists():
                media_files = sum(1 for p in media_dir.rglob("*") if p.is_file())
            files_dir = run_dir / "files"
            if files_dir.exists():
                save_links = sum(1 for p in files_dir.rglob("*") if p.is_symlink())

            total_s = t_total_end - t_total_start

            result = {
                "skip_store": skip_store,
                "steps": steps,
                "keys": keys,
                "scalar_records": steps * keys,
                "media_records": num_image,
                "save_files": save_count,
                "producer_s": round(t_producer_end - t_total_start, 4),
                "total_s": round(total_s, 4),
                "finish_ms": round(finish_ms, 2),
                "producer_rec_per_s": round(steps * keys / (t_producer_end - t_total_start), 1),
                "e2e_rec_per_s": round(steps * keys / total_s, 1),
                "run_file_bytes": run_file_bytes,
                "media_files": media_files,
                "save_links": save_links,
                "log": _latency_stats(log_lat),
            }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    assert result, "worker produced no result"
    print(RESULT_MARK + json.dumps(result), flush=True)


def spawn_case(skip_store: bool, steps: int, keys: int, save_count: int) -> Dict[str, Any]:
    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--worker",
        "--skip-store" if skip_store else "--persist",
        "--steps",
        str(steps),
        "--keys",
        str(keys),
        "--save-count",
        str(save_count),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    for line in proc.stdout.splitlines():
        if line.startswith(RESULT_MARK):
            return json.loads(line[len(RESULT_MARK) :])
    sys.stderr.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    raise RuntimeError(f"Worker failed: skip_store={skip_store}, rc={proc.returncode}")


# ===========================================================================
# 报告
# ===========================================================================


def best_of(skip_store: bool, steps: int, keys: int, save_count: int, repeats: int) -> Dict[str, Any]:
    """重复运行并返回 total_s 最优的一次（降低调度噪声）。"""
    best: Optional[Dict[str, Any]] = None
    for _ in range(repeats):
        r = spawn_case(skip_store=skip_store, steps=steps, keys=keys, save_count=save_count)
        if best is None or r["total_s"] < best["total_s"]:
            best = r
    assert best is not None
    return best


def _print_report(persist: Dict[str, Any], skip: Dict[str, Any]) -> None:
    print("\n" + "=" * 92)
    print("  SwanLab Benchmark: core.skip_store end-to-end (online + mocked HTTP)")
    print("=" * 92)
    print(
        f"  workload: {persist['steps']} steps x {persist['keys']} keys "
        f"= {persist['scalar_records']:,} scalars, {persist['media_records']} media, "
        f"{persist['save_files']} saves"
    )
    print("-" * 92)
    print(f"  {'metric':<28}{'persist':>18}{'skip_store':>18}{'delta':>20}")
    print("-" * 92)

    def row(label: str, a: Any, b: Any, lower_better: bool = True) -> None:
        if isinstance(a, float) and isinstance(b, float) and a != 0:
            pct = (a - b) / a * 100.0
            delta = f"{pct:+.1f}%"
        else:
            delta = "-"
        print(f"  {label:<28}{str(a):>18}{str(b):>18}{delta:>20}")

    row("producer (s)", persist["producer_s"], skip["producer_s"])
    row("total wall (s)", persist["total_s"], skip["total_s"])
    row("finish drain (ms)", persist["finish_ms"], skip["finish_ms"])
    row("producer throughput (rec/s)", persist["producer_rec_per_s"], skip["producer_rec_per_s"])
    row("e2e throughput (rec/s)", persist["e2e_rec_per_s"], skip["e2e_rec_per_s"])
    print("-" * 92)
    print("  run.log call (us)  [主线程耗时；标量为主，不含后台上传]")
    for stat in ("mean", "p50", "p95", "p99"):
        print(f"  {'  ' + stat:<28}{persist['log'][stat]:>18}{skip['log'][stat]:>18}")
    print("-" * 92)
    row("run-*.swanlab bytes", persist["run_file_bytes"], skip["run_file_bytes"])
    row("media files", persist["media_files"], skip["media_files"])
    row("save links", persist["save_links"], skip["save_links"])
    print("-" * 92)
    print("  note: producer 阶段含 BackgroundConsumer 的本地 store 工作。skip_store 把它从")
    print("        磁盘 I/O 变成内存操作，主线程与消费者线程的 GIL 竞争可能反而拉长 producer")
    print("        墙钟；真正的收益体现在 finish 排空耗时与端到端总耗时上。")
    print("=" * 92)


# ===========================================================================
# Pytest 入口
# ===========================================================================


def test_bench_skip_store_e2e():
    """online + mock HTTP 下端到端对比 skip_store 的吞吐量与延时。"""
    persist = best_of(skip_store=False, steps=STEPS, keys=KEYS, save_count=SAVE_COUNT, repeats=REPEATS)
    skip = best_of(skip_store=True, steps=STEPS, keys=KEYS, save_count=SAVE_COUNT, repeats=REPEATS)

    _print_report(persist, skip)

    # 正确性：persist 有本地产物，skip 完全没有
    assert persist["run_file_bytes"] > 0
    assert persist["media_files"] > 0
    assert persist["save_links"] > 0
    assert skip["run_file_bytes"] == 0
    assert skip["media_files"] == 0
    assert skip["save_links"] == 0
    # 场景一致性
    assert persist["scalar_records"] == skip["scalar_records"]
    assert persist["media_records"] == skip["media_records"]


# ===========================================================================
# CLI
# ===========================================================================


def main() -> None:
    parser = argparse.ArgumentParser(description="SwanLab: core.skip_store end-to-end benchmark")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--skip-store", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--persist", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--steps", type=int, default=STEPS)
    parser.add_argument("--keys", type=int, default=KEYS)
    parser.add_argument("--save-count", type=int, default=SAVE_COUNT)
    parser.add_argument("--repeats", type=int, default=REPEATS)
    args = parser.parse_args()

    if args.worker:
        run_worker(skip_store=args.skip_store, steps=args.steps, keys=args.keys, save_count=args.save_count)
        return

    persist = best_of(
        skip_store=False, steps=args.steps, keys=args.keys, save_count=args.save_count, repeats=args.repeats
    )
    skip = best_of(skip_store=True, steps=args.steps, keys=args.keys, save_count=args.save_count, repeats=args.repeats)
    _print_report(persist, skip)


if __name__ == "__main__":
    main()
