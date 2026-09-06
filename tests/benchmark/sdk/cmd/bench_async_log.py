"""
@author: cunyue, nexisato
@file: bench_async_log.py
@time: 2026/9/4
@description: swanlab.log vs swanlab.async_log (threading 模式) 性能对比基准测试

对比维度：
  1. 上报模式：
     - no_sdk:    纯计算基线（无任何 SwanLab 开销，仅 macro 作为 0% 黄金基准）
     - sync_log:  swanlab.log(data, step=step) 同步上报
     - async_log: swanlab.async_log(lambda: data, step=step, mode="threading") 异步上报
  2. 指标规模：
     - 10 keys:   常规训练常用标量（loss, acc, lr, 层权重/梯度模长等）
     - 100 keys:  宽/层级命名空间标量（大模型各层细节监控）
  3. 评测类型：
     - micro: 纯主线程调用耗时（Mean / P50 / P95 / P99 延迟，Throughput QPS）
     - macro: 模拟训练循环（Step 耗时分布，相对 No-SDK 的 Overhead % 开销占比）

运行方式：
  # 1. 作为 pytest 基准测试运行：
  uv run pytest tests/benchmark/sdk/cmd/bench_async_log.py -v -s

  # 2. 作为独立脚本运行（支持完整 CLI 参数）：
  uv run python tests/benchmark/sdk/cmd/bench_async_log.py
  uv run python tests/benchmark/sdk/cmd/bench_async_log.py --steps 500 --workload-ms 10
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
from typing import Any, Callable, Dict, List, Optional

RESULT_MARK = "##BENCH_JSON## "


# ===========================================================================
# 指标数据生成器 (10 keys vs 100 keys)
# ===========================================================================


def generate_payload(num_keys: int, step: int) -> Dict[str, float]:
    """生成指定 key 数量的指标字典，包含层级路径。"""
    if num_keys == 10:
        return {
            "train/loss": 0.5 + 0.001 * (step % 100),
            "train/accuracy": 0.85 + 0.0005 * (step % 100),
            "train/learning_rate": 1e-4,
            "train/grad_norm": 1.25,
            "train/epoch": step / 100.0,
            "layer1/weight_norm": 0.35,
            "layer1/bias_norm": 0.05,
            "layer2/weight_norm": 0.45,
            "layer2/bias_norm": 0.06,
            "val/loss": 0.52 + 0.001 * (step % 100),
        }
    elif num_keys == 100:
        data: Dict[str, float] = {}
        data["train/loss"] = 0.5 + 0.001 * (step % 100)
        data["train/accuracy"] = 0.85 + 0.0005 * (step % 100)
        data["train/lr"] = 1e-4
        data["train/epoch"] = step / 100.0
        # 其余 96 个指标按 12 个模块 × 8 个统计量组织
        stats = ["w_mean", "w_std", "b_mean", "b_std", "grad_norm", "grad_std", "act_mean", "act_std"]
        for module_idx in range(12):
            for stat in stats:
                data[f"layers/block_{module_idx:02d}/{stat}"] = 0.01 * (step % 50) + module_idx * 0.1
        return data
    else:
        return {f"metric_{i:03d}": 1.0 + i * 0.01 for i in range(num_keys)}


# ===========================================================================
# 统计辅助函数
# ===========================================================================


def percentile(sorted_data: List[float], pct: float) -> float:
    """计算百分位数，pct ∈ [0, 100]"""
    if not sorted_data:
        return 0.0
    idx = int(len(sorted_data) * pct / 100.0)
    return sorted_data[min(idx, len(sorted_data) - 1)]


def simulate_training_workload(workload_ms: float) -> None:
    """模拟真实训练计算（释放 GIL 的等待 + 少量 CPU 计算）。"""
    if workload_ms <= 0:
        return
    # 模拟真实混合计算：80% sleep (模拟 GPU forward/backward C 扩展), 20% CPU 占用
    sleep_s = (workload_ms * 0.8) / 1000.0
    time.sleep(sleep_s)
    cpu_deadline = time.perf_counter() + (workload_ms * 0.2) / 1000.0
    val = 1.0001
    while time.perf_counter() < cpu_deadline:
        val = (val * 1.000001) % 1000.0


# ===========================================================================
# 子进程 Worker 实现（彻底隔离环境，消除全局单例和线程残留污染）
# ===========================================================================


def run_worker(
    mode: str,
    num_keys: int,
    test_type: str,
    steps: int,
    warmup: int,
    workload_ms: float,
) -> None:
    """在独立子进程中执行单次压测。

    :param mode: "no_sdk" | "sync_log" | "async_log"
    :param num_keys: 10 | 100
    :param test_type: "micro" | "macro"
    :param steps: 测量步骤数
    :param warmup: 预热步骤数（耗时不计入统计）
    :param workload_ms: macro 模式下的每步模拟训练耗时 (ms)
    """
    tmp_log_dir = tempfile.mkdtemp(prefix="swanlab_bench_")
    call_latencies_ns: List[int] = []
    step_latencies_ns: List[int] = []

    # 上报调用闭包，在 init 分支内绑定到具体实现，热循环内不做模式分支
    log_fn: Optional[Callable[[Dict[str, float], int], Any]] = None
    finish_fn: Optional[Callable[[], Any]] = None

    try:
        if mode != "no_sdk":
            import swanlab

            # 统一使用 local 模式排除网络波动，设置 log_level 减少控制台干扰
            swanlab.init(
                project="benchmarks",
                experiment_name=f"bench_{mode}_{num_keys}k_{test_type}",
                mode="local",
                log_level="warning",
                log_dir=tmp_log_dir,
            )
            finish_fn = swanlab.finish
            if mode == "sync_log":

                def sync_log_fn(d: Dict[str, float], s: int) -> None:
                    swanlab.log(d, step=s)

                log_fn = sync_log_fn
            else:

                def async_log_fn(d: Dict[str, float], s: int) -> Any:
                    return swanlab.async_log(lambda d=d: d, step=s, mode="threading")

                log_fn = async_log_fn

        # ----------------------------------
        # 预热阶段 (Warmup)
        # ----------------------------------
        for w in range(warmup):
            data = generate_payload(num_keys, w)
            if test_type == "macro":
                simulate_training_workload(workload_ms)
            if log_fn is not None:
                log_fn(data, w)

        # ----------------------------------
        # 测量阶段 (Measurement)
        # ----------------------------------
        t_total_start = time.perf_counter_ns()

        for s in range(warmup, warmup + steps):
            data = generate_payload(num_keys, s)

            t_step_start = time.perf_counter_ns()

            if test_type == "macro":
                simulate_training_workload(workload_ms)

            t_call_start = time.perf_counter_ns()
            if log_fn is not None:
                log_fn(data, s)
            t_call_end = time.perf_counter_ns()

            t_step_end = time.perf_counter_ns()

            call_latencies_ns.append(t_call_end - t_call_start)
            step_latencies_ns.append(t_step_end - t_step_start)

        # 主循环耗时到此为止；finish 的队列排空耗时单独统计，不计入吞吐量
        t_total_end = time.perf_counter_ns()

        t_finish_start = time.perf_counter_ns()
        if finish_fn is not None:
            finish_fn()
        t_finish_end = time.perf_counter_ns()

        # ----------------------------------
        # 统计分析
        # ----------------------------------
        call_latencies_us = [v / 1e3 for v in call_latencies_ns]
        step_latencies_ms = [v / 1e6 for v in step_latencies_ns]

        call_latencies_us_sorted = sorted(call_latencies_us)
        step_latencies_ms_sorted = sorted(step_latencies_ms)

        total_time_s = (t_total_end - t_total_start) / 1e9
        finish_time_ms = (t_finish_end - t_finish_start) / 1e6

        result = {
            "mode": mode,
            "num_keys": num_keys,
            "test_type": test_type,
            "steps": steps,
            "warmup": warmup,
            "workload_ms": workload_ms,
            "total_time_s": round(total_time_s, 4),
            "finish_time_ms": round(finish_time_ms, 2),
            "throughput_qps": round(steps / total_time_s, 1),
            # 单次调用延迟 (μs)
            "call_mean_us": round(sum(call_latencies_us) / len(call_latencies_us), 2),
            "call_p50_us": round(percentile(call_latencies_us_sorted, 50), 2),
            "call_p95_us": round(percentile(call_latencies_us_sorted, 95), 2),
            "call_p99_us": round(percentile(call_latencies_us_sorted, 99), 2),
            # Step 耗时 (ms)
            "step_mean_ms": round(sum(step_latencies_ms) / len(step_latencies_ms), 3),
            "step_p50_ms": round(percentile(step_latencies_ms_sorted, 50), 3),
            "step_p95_ms": round(percentile(step_latencies_ms_sorted, 95), 3),
            "step_p99_ms": round(percentile(step_latencies_ms_sorted, 99), 3),
        }
        print(RESULT_MARK + json.dumps(result), flush=True)

    finally:
        shutil.rmtree(tmp_log_dir, ignore_errors=True)


# ===========================================================================
# 父进程调用与结果汇总
# ===========================================================================


def spawn_case(
    mode: str,
    num_keys: int,
    test_type: str,
    steps: int,
    warmup: int,
    workload_ms: float,
) -> Dict[str, Any]:
    """在子进程中运行单个测试用例并解析结果。"""
    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--worker",
        "--mode",
        mode,
        "--keys",
        str(num_keys),
        "--test-type",
        test_type,
        "--steps",
        str(steps),
        "--warmup",
        str(warmup),
        "--workload-ms",
        str(workload_ms),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    payload: Optional[Dict[str, Any]] = None
    for line in proc.stdout.splitlines():
        if line.startswith(RESULT_MARK):
            payload = json.loads(line[len(RESULT_MARK) :])
            break

    if payload is None:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Worker failed: mode={mode}, keys={num_keys}, rc={proc.returncode}")

    return payload


def print_report(
    micro_results: List[Dict[str, Any]],
    macro_results: List[Dict[str, Any]],
) -> None:
    """打印格式化的高亮性能对比报告。"""
    print("\n" + "=" * 108)
    print("  SwanLab Benchmark: swanlab.log vs swanlab.async_log (threading 模式)")
    print("=" * 108)

    # 1. Micro Benchmark: 纯主线程单次调用延迟对比
    print("\n[Part 1] Micro-Benchmark: 纯主线程调用耗时与吞吐量 (无模拟训练开销)")
    print("-" * 108)
    print(
        f"  {'Keys':<8} {'Mode':<14} {'Call Mean(μs)':>14} {'P50(μs)':>10} {'P95(μs)':>10} {'P99(μs)':>10} {'Throughput(QPS)':>16} {'Finish(ms)':>12}"
    )
    print("-" * 108)

    for r in micro_results:
        print(
            f"  {r['num_keys']:<8} "
            f"{r['mode']:<14} "
            f"{r['call_mean_us']:>14.1f} "
            f"{r['call_p50_us']:>10.1f} "
            f"{r['call_p95_us']:>10.1f} "
            f"{r['call_p99_us']:>10.1f} "
            f"{r['throughput_qps']:>16.1f} "
            f"{r['finish_time_ms']:>12.1f}"
        )
    print("-" * 108)

    # 2. Macro Benchmark: 模拟训练循环耗时与开销占比对比
    print("\n[Part 2] Macro-Benchmark: 模拟训练 Step 耗时分布与开销占比 (%)")
    print("-" * 108)
    print(
        f"  {'Keys':<8} {'Mode':<14} {'Step Mean(ms)':>14} {'P95(ms)':>10} {'P99(ms)':>10} {'Overhead(%)':>13} {'vs sync_log':>16} {'Finish(ms)':>12}"
    )
    print("-" * 108)

    # 按 num_keys 分组计算开销百分比
    grouped_keys: Dict[int, Dict[str, Dict[str, Any]]] = {}
    for r in macro_results:
        k = r["num_keys"]
        grouped_keys.setdefault(k, {})[r["mode"]] = r

    for k, group in grouped_keys.items():
        base = group.get("no_sdk")
        base_ms = base["step_mean_ms"] if base else 1.0
        sync_r = group.get("sync_log")
        sync_overhead = (sync_r["step_mean_ms"] - base_ms) / base_ms * 100.0 if sync_r and base else 0.0

        for mode_name in ("no_sdk", "sync_log", "async_log"):
            if mode_name not in group:
                continue
            r = group[mode_name]
            cur_ms = r["step_mean_ms"]
            overhead_pct = (cur_ms - base_ms) / base_ms * 100.0 if base else 0.0

            if mode_name == "no_sdk":
                overhead_str = "0.00% (基准)"
                vs_sync_str = "-"
            elif mode_name == "sync_log":
                overhead_str = f"+{overhead_pct:.2f}%"
                vs_sync_str = "基准 (100%)"
            else:
                overhead_str = f"+{overhead_pct:.2f}%"
                if sync_overhead > 0.001:
                    reduction = (sync_overhead - overhead_pct) / sync_overhead * 100.0
                    vs_sync_str = f"开销减少 {reduction:.1f}%"
                else:
                    vs_sync_str = "相当"

            print(
                f"  {k:<8} "
                f"{mode_name:<14} "
                f"{r['step_mean_ms']:>14.3f} "
                f"{r['step_p95_ms']:>10.3f} "
                f"{r['step_p99_ms']:>10.3f} "
                f"{overhead_str:>13} "
                f"{vs_sync_str:>16} "
                f"{r['finish_time_ms']:>12.1f}"
            )
        print("-" * 108)

    print("=" * 108 + "\n")


def run_full_benchmark(
    micro_steps: int = 1000,
    macro_steps: int = 300,
    workload_ms: float = 10.0,
    warmup: int = 50,
) -> Dict[str, Any]:
    """执行完整的 10-key 和 100-key 对比测试矩阵。"""
    keys_list = [10, 100]
    modes_micro = ["sync_log", "async_log"]
    modes_macro = ["no_sdk", "sync_log", "async_log"]

    micro_results: List[Dict[str, Any]] = []
    macro_results: List[Dict[str, Any]] = []

    # 1. 运行 Micro-benchmark
    for k in keys_list:
        for m in modes_micro:
            sys.stderr.write(f"[Bench] Micro test: keys={k}, mode={m}, steps={micro_steps}...\n")
            res = spawn_case(m, k, "micro", micro_steps, warmup, workload_ms=0.0)
            micro_results.append(res)

    # 2. 运行 Macro-benchmark
    for k in keys_list:
        for m in modes_macro:
            sys.stderr.write(
                f"[Bench] Macro test: keys={k}, mode={m}, steps={macro_steps}, workload={workload_ms}ms...\n"
            )
            res = spawn_case(m, k, "macro", macro_steps, warmup, workload_ms=workload_ms)
            macro_results.append(res)

    print_report(micro_results, macro_results)
    return {"micro": micro_results, "macro": macro_results}


# ===========================================================================
# Pytest 测试入口
# ===========================================================================


def test_bench_log_vs_async_log():
    """Pytest 基准测试套件集成入口。"""
    # pytest 下适当缩短步数以保证测试快速通过，同时保留统计精度
    results = run_full_benchmark(
        micro_steps=500,
        macro_steps=200,
        workload_ms=5.0,
        warmup=30,
    )
    # 验证测试完整性
    assert len(results["micro"]) == 4  # 2 keys * 2 modes
    assert len(results["macro"]) == 6  # 2 keys * 3 modes

    # 验证 100-key 场景下 async_log 主线程单次调用延迟显著优于 sync_log
    micro_100 = {r["mode"]: r for r in results["micro"] if r["num_keys"] == 100}
    assert micro_100["async_log"]["call_mean_us"] < micro_100["sync_log"]["call_mean_us"], (
        "async_log should have lower main-thread latency than sync_log on 100 keys"
    )


# ===========================================================================
# CLI 独立运行入口
# ===========================================================================


def main():
    parser = argparse.ArgumentParser(description="SwanLab: swanlab.log vs swanlab.async_log Benchmark")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--mode", type=str, choices=["no_sdk", "sync_log", "async_log"], default="sync_log")
    parser.add_argument("--keys", type=int, default=10, help="指标 key 数量 (10 或 100)")
    parser.add_argument("--test-type", type=str, choices=["micro", "macro"], default="micro")
    parser.add_argument("--steps", type=int, default=1000, help="测试 step 数")
    parser.add_argument("--warmup", type=int, default=50, help="预热 step 数")
    parser.add_argument("--workload-ms", type=float, default=10.0, help="模拟训练每步耗时 (ms)")
    parser.add_argument("--out", type=str, default=None, help="导出结果 JSON 路径")
    args = parser.parse_args()

    if args.worker:
        run_worker(
            mode=args.mode,
            num_keys=args.keys,
            test_type=args.test_type,
            steps=args.steps,
            warmup=args.warmup,
            workload_ms=args.workload_ms,
        )
        return

    results = run_full_benchmark(
        micro_steps=args.steps,
        macro_steps=max(200, args.steps // 3),
        workload_ms=args.workload_ms,
        warmup=args.warmup,
    )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
        print(f"Results saved to {args.out}")


if __name__ == "__main__":
    main()
