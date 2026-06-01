"""
@author: nexisato
@file: bench_store_fsync.py
@time: 2026/6/1
@description: DataStoreWriter fsync 优化前后性能对比基准测试

直接调用 DataStoreWriter.write() 写入固定大小 payload，
隔离 SDK 运行时开销，精确测量磁盘 IO 差异。

关注指标：
  1. 总耗时
  2. 吞吐量 (rec/s)
  3. 单次 write 延迟 (mean / p50 / p95 / p99 / max)
  4. 文件大小（验证一致性）
  5. 回读记录数（验证数据完整性）

用法：
  uv run pytest tests/benchmark/sdk/internal/core_python/store/bench_store_fsync.py -v -s
"""

import time
from typing import List

from swanlab.sdk.internal.core_python.store import DataStoreReader, DataStoreWriter

# ===========================================================================
# 参数配置
# ===========================================================================

NUM_RECORDS = 100_000
PAYLOAD_SIZE = 120  # 模拟标量 protobuf 序列化后大小


# ===========================================================================
# 辅助函数
# ===========================================================================


def _percentile(sorted_data: List[float], pct: float) -> float:
    """计算百分位，pct ∈ [0, 100]"""
    if not sorted_data:
        return 0.0
    idx = int(len(sorted_data) * pct / 100)
    return sorted_data[min(idx, len(sorted_data) - 1)]


# ===========================================================================
# Benchmark
# ===========================================================================


def test_bench_raw_store(tmp_path):
    """
    DataStoreWriter 裸写性能基准测试。

    写入 10,000 条 120 bytes 的 record，测量延迟分布和吞吐量，
    最后回读验证数据完整性。
    """
    store_path = tmp_path / "raw_bench.swanlab"

    # 预生成 protobuf-like payload（带一些变化）
    payloads = [bytes([i % 256]) * PAYLOAD_SIZE for i in range(min(NUM_RECORDS, 256))]
    if NUM_RECORDS > 256:
        payloads = payloads * (NUM_RECORDS // 256 + 1)
    payloads = payloads[:NUM_RECORDS]

    # ---- 写入 benchmark ----
    writer = DataStoreWriter()
    writer.open(str(store_path))

    latencies: List[float] = []
    t_total_start = time.perf_counter()

    for i in range(NUM_RECORDS):
        t0 = time.perf_counter()
        writer.write(payloads[i])
        latencies.append(time.perf_counter() - t0)

    writer.close()
    t_total_end = time.perf_counter()

    # ---- 统计 ----
    latencies.sort()
    total_time = t_total_end - t_total_start
    mean_latency_us = sum(latencies) / len(latencies) * 1e6

    file_size = store_path.stat().st_size

    # 回读验证
    reader = DataStoreReader()
    reader.open(str(store_path))
    read_count = sum(1 for _ in reader)
    reader.close()

    result = {
        "num_records": NUM_RECORDS,
        "payload_size_bytes": PAYLOAD_SIZE,
        "total_time_sec": round(total_time, 4),
        "throughput_rec_per_sec": round(NUM_RECORDS / total_time, 1),
        "write_mean_us": round(mean_latency_us, 2),
        "write_p50_us": round(_percentile(latencies, 50) * 1e6, 2),
        "write_p95_us": round(_percentile(latencies, 95) * 1e6, 2),
        "write_p99_us": round(_percentile(latencies, 99) * 1e6, 2),
        "write_max_us": round(latencies[-1] * 1e6, 2),
        "file_size_bytes": file_size,
        "reader_verify_count": read_count,
    }

    print("\n" + "=" * 60)
    print("  Raw DataStoreWriter Benchmark")
    print("=" * 60)
    for k, v in result.items():
        print(f"  {k:30s}: {v}")
    print("=" * 60)

    # 回读数应 == 写入数
    assert read_count == NUM_RECORDS, f"Data integrity check failed: wrote {NUM_RECORDS}, read {read_count}"
