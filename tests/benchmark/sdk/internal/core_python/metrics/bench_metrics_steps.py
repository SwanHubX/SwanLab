"""
@author: nexisato
@file: bench_metrics_steps.py
@time: 2026/8/18
@description: BaseMetric.steps 去重容器优化前后性能对比基准测试

对比 Set[int]（旧实现）与 pyroaring BitMap64 + 周期 run_optimize（现行实现）
在指定 key 数 × step 数下的内存占用与吞吐量。

为保证 set 基线为旧版生产行为的逐字复刻、bitmap 走现行真实路径（RunMetrics
→ ScalarMetric.try_accept_step），每个 (实现, 工作负载) 组合在独立子进程中
运行，通过 RSS 增量测量常驻内存，避免同进程交叉污染。

工作负载：
  1. sequential: 每个 key 顺序写入 1..NUM_STEPS（真实训练形态）
  2. backfill:   顺序写入 + 周期性重复/回填探测（验证去重拒绝路径）
  3. sparse:     每个 key 在 [1, 2^48) 内随机不重复 step（极限bad case下的病态稀疏）

关注指标：
  1. 吞吐量 (ops/s)
  2. 常驻内存增量 (RSS delta) 与每 key 内存
  3. 序列化大小（bitmap: run 容器收敛验证）
  4. accepted/rejected 计数与 set 基线差分一致（正确性）

用法：
  uv run pytest tests/benchmark/sdk/internal/core_python/metrics/bench_metrics_steps.py -v -s
"""

import json
import os
import random
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

# ===========================================================================
# 参数配置
# ===========================================================================

NUM_KEYS = 200
NUM_STEPS = 5_000
SPARSE_RANGE = 1 << 48  # 随机 step 上界（开放，bit 空间上限 _MAX_STEP）

IMPLS = ("set", "bitmap")
WORKLOADS = ("sequential", "backfill", "sparse")


# ===========================================================================
# 辅助函数
# ===========================================================================


def _rss_kib() -> int:
    """当前进程常驻内存（KiB），macOS / Linux 通用"""
    out = subprocess.check_output(["ps", "-o", "rss=", "-p", str(os.getpid())], text=True)
    return int(out.strip())


def _gen_ops(workload: str, key_index: int, num_steps: int) -> List[int]:
    """生成单个 key 的 step 操作序列（含去重/回填探测）"""
    if workload == "sequential":
        return list(range(1, num_steps + 1))
    if workload == "backfill":
        ops: List[int] = []
        for s in range(1, num_steps + 1):
            ops.append(s)
            if s % 10 == 0:
                ops.append(s // 2)  # 重复：必拒绝
                if s > 5:
                    ops.append(s - 5)  # 回填：必拒绝
        return ops
    if workload == "sparse":
        rng = random.Random(f"{workload}:{key_index}")
        return rng.sample(range(1, SPARSE_RANGE), num_steps)
    raise ValueError(f"unknown workload: {workload}")


class _SetStepTracker:
    """旧实现逐字复刻（git 历史中 BaseMetric.try_accept_step 的 set 版本）"""

    __slots__ = ("steps", "min_step")

    def __init__(self) -> None:
        self.steps: set = set()
        self.min_step: int = -1

    def try_accept_step(self, step: int) -> bool:
        if step <= self.min_step:
            return False
        if step in self.steps:
            return False
        self.steps.add(step)
        return True


# ===========================================================================
# 子进程模式：单个 (impl, workload) 组合的测量
# ===========================================================================


def _child_main(impl: str, workload: str) -> None:
    from swanlab.sdk.internal.core_python.metrics import RunMetrics
    from swanlab.sdk.internal.core_python.pkg import builder

    def make_trackers() -> Tuple[List[Any], Optional[int]]:
        """返回 (tracker 列表, 序列化字节样本所在下标或 None)"""
        if impl == "set":
            return [_SetStepTracker() for _ in range(NUM_KEYS)], None
        run_metrics = RunMetrics()
        return (
            [
                run_metrics.define_scalar(key=f"k{i}", column=builder.build_resume_column(f"k{i}"))
                for i in range(NUM_KEYS)
            ],
            0,
        )

    # ---- 预热：导入 protobuf / 编译热点 / 触发模块级懒加载 ----
    warm = make_trackers()
    for t in warm[0][:1]:
        for s in range(1, 257):
            t.try_accept_step(s)
    del warm

    import gc

    gc.collect()
    rss0 = _rss_kib()

    # ---- 正式测量 ----
    trackers, ser_index = make_trackers()
    accepted = rejected = 0
    t_total_start = time.perf_counter()

    for key_index, tracker in enumerate(trackers):
        ops = _gen_ops(workload, key_index, NUM_STEPS)
        for s in ops:
            if tracker.try_accept_step(s):
                accepted += 1
            else:
                rejected += 1

    t_total_end = time.perf_counter()

    gc.collect()
    rss1 = _rss_kib()

    total_ops = accepted + rejected
    total_time = t_total_end - t_total_start

    serialize_bytes: Optional[int] = None
    if ser_index is not None:
        serialize_bytes = len(trackers[ser_index].steps.serialize())

    stats = {
        "impl": impl,
        "workload": workload,
        "num_keys": NUM_KEYS,
        "num_steps": NUM_STEPS,
        "ops": total_ops,
        "accepted": accepted,
        "rejected": rejected,
        "total_time_sec": round(total_time, 4),
        "ops_per_sec": round(total_ops / total_time, 1),
        "rss_delta_mb": round((rss1 - rss0) / 1024, 2),
        "per_key_kb": round((rss1 - rss0) / NUM_KEYS, 2),
        "serialize_key0_bytes": serialize_bytes,
    }
    print(json.dumps(stats))


# ===========================================================================
# Benchmark
# ===========================================================================


_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../.."))


def _run_case(impl: str, workload: str) -> Dict[str, Any]:
    """独立子进程运行单个组合，返回解析后的统计"""
    proc = subprocess.run(
        [sys.executable, os.path.abspath(__file__), impl, workload],
        capture_output=True,
        text=True,
        cwd=_REPO_ROOT,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"child {impl}/{workload} failed:\n{proc.stderr[-2000:]}")
    json_lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if not json_lines:
        raise RuntimeError(f"child {impl}/{workload} produced no JSON:\n{proc.stdout}\n{proc.stderr[-2000:]}")
    return json.loads(json_lines[-1])


def _print_pair(workload: str, r_set: Dict[str, Any], r_bm: Dict[str, Any]) -> None:
    mem_ratio = r_set["per_key_kb"] / r_bm["per_key_kb"] if r_bm["per_key_kb"] > 0 else float("inf")
    tput_ratio = r_set["ops_per_sec"] / r_bm["ops_per_sec"] if r_bm["ops_per_sec"] > 0 else float("inf")
    rows = [
        ("ops", f"{r_set['ops']}", f"{r_bm['ops']}"),
        (
            "accepted / rejected",
            f"{r_set['accepted']} / {r_set['rejected']}",
            f"{r_bm['accepted']} / {r_bm['rejected']}",
        ),
        ("total_time_sec", f"{r_set['total_time_sec']}", f"{r_bm['total_time_sec']}"),
        ("ops_per_sec", f"{r_set['ops_per_sec']}", f"{r_bm['ops_per_sec']}"),
        ("rss_delta_mb", f"{r_set['rss_delta_mb']}", f"{r_bm['rss_delta_mb']}"),
        ("per_key_kb", f"{r_set['per_key_kb']}", f"{r_bm['per_key_kb']}"),
        ("serialize_key0_bytes", "-", f"{r_bm['serialize_key0_bytes']}"),
        ("mem_ratio (set/bitmap)", "-", f"{mem_ratio:.3g}x"),
        ("tput_ratio (set/bitmap)", "-", f"{tput_ratio:.3g}x"),
    ]
    print("\n" + "=" * 74)
    print(f"  Steps Tracker Benchmark: {workload}  (keys={NUM_KEYS}, steps={NUM_STEPS})")
    print("=" * 74)
    print(f"  {'metric':28s}: {'set':>16s} | {'bitmap':>16s}")
    print("-" * 74)
    for name, v_set, v_bm in rows:
        print(f"  {name:28s}: {v_set:>16s} | {v_bm:>16s}")
    print("=" * 74)


def test_bench_steps_tracker():
    """
    Set[int] vs BitMap64 step 去重容器对比基准。

    每个 (实现, 工作负载) 在独立子进程中运行，测量 RSS 增量与吞吐量，
    最后验证 bitmap 与 set 基线的 accepted/rejected 差分一致。
    """
    results = {}
    for workload in WORKLOADS:
        for impl in IMPLS:
            results[(workload, impl)] = _run_case(impl, workload)

    for workload in WORKLOADS:
        r_set, r_bm = results[(workload, "set")], results[(workload, "bitmap")]
        # 差分一致性：两种实现对同一操作序列必须产生相同的接受/拒绝判定
        assert r_bm["accepted"] == r_set["accepted"], (
            f"[{workload}] accepted mismatch: set={r_set['accepted']}, bitmap={r_bm['accepted']}"
        )
        assert r_bm["rejected"] == r_set["rejected"], (
            f"[{workload}] rejected mismatch: set={r_set['rejected']}, bitmap={r_bm['rejected']}"
        )
        _print_pair(workload, r_set, r_bm)

    # sequential 全部接受，backfill 拒绝数为接受数的 0.2 倍
    seq = results[("sequential", "bitmap")]
    assert seq["accepted"] == NUM_KEYS * NUM_STEPS
    bf = results[("backfill", "bitmap")]
    assert bf["accepted"] == NUM_KEYS * NUM_STEPS
    assert bf["rejected"] == NUM_KEYS * (NUM_STEPS // 10) * 2


if __name__ == "__main__":
    _child_main(sys.argv[1], sys.argv[2])
