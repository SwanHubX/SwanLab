"""
@author: nexisato
@file: bench_store_skip.py
@time: 2026/9/11
@description: core.skip_store 本地持久化开销基准测试

对比 online mode 落盘（skip_store=False）与完全跳过持久化（skip_store=True）
在相同写入负载下的性能，量化“跳过 protobuf 序列化 + DataStore 落盘”带来的收益。

两个场景均走生产路径 ``CorePython._store_records``：
  - skip_store=False: 逐条 ``Record.SerializeToString()`` 后写入 LevelDB log
  - skip_store=True : 仅累加未持久化计数，不序列化、不创建文件

关注指标：
  1. 总耗时 (s)
  2. 吞吐量 (rec/s)
  3. 单条延迟 (us)
  4. 加速比 (normal / skip)
  5. 数据完整性（normal 回读记录数 == 写入数；skip 不产生文件）

用法：
  uv run pytest tests/benchmark/sdk/internal/core_python/store/bench_store_skip.py -v -s
"""

import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord, ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.core import CorePython
from swanlab.sdk.internal.core_python.store import DataStoreReader, DataStoreWriter

# ===========================================================================
# 参数配置
# ===========================================================================

# 与 bench_metrics_steps.py 对齐：200 key × 5000 step，共 1,000,000 条记录
# （总记录数同时与 bench_store_fsync.py 的 NUM_RECORDS = 1,000,000 一致）
NUM_KEYS = 200
NUM_STEPS = 5_000
NUM_RECORDS = NUM_KEYS * NUM_STEPS
BATCH_SIZE = 100  # 与 BackgroundConsumer 的 scalar 批大小一致
REPEATS = 3  # 每场景重复次数，取最优以降低噪声


# ===========================================================================
# 辅助函数
# ===========================================================================


def _make_records(num_keys: int, num_steps: int) -> List[Record]:
    """预生成 step × key 的标量 Record，避免把构建开销算进持久化基准。

    key 命名与 bench_metrics_steps.py 的 ``k{i}`` 一致，便于横向对比。
    """
    records: List[Record] = []
    num = 0
    for step in range(1, num_steps + 1):
        for k in range(num_keys):
            num += 1
            records.append(
                Record(
                    num=num,
                    scalar=ScalarRecord(
                        key=f"k{k}",
                        step=step,
                        type=ColumnType.COLUMN_TYPE_SCALAR,
                        value=ScalarValue(number=0.5 + (step % 1000) * 1e-4),
                    ),
                )
            )
    return records


def _make_core(run_dir: Path, skip_store: bool) -> CorePython:
    """构造仅用于持久化路径的 CorePython（不启动 transport / heartbeat）。"""
    core = CorePython(mode="online")
    core._ctx = CoreContext(
        config=CoreConfig(
            run_id="bench-run",
            run_dir=run_dir,
            section_rule=0,
            record_batch=BATCH_SIZE,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
            skip_store=skip_store,
        )
    )
    core._store = DataStoreWriter(skip=skip_store)
    core._store.open(str(core._ctx.run_file))
    return core


def _write_all(core: CorePython, records: List[Record]) -> float:
    """按 BATCH_SIZE 分批调用生产路径 _store_records，返回耗时（秒）。"""
    start = time.perf_counter()
    for offset in range(0, len(records), BATCH_SIZE):
        core._store_records(records[offset : offset + BATCH_SIZE])
    return time.perf_counter() - start


def _run_scenario(root: Path, skip_store: bool, records: List[Record]) -> Tuple[float, CorePython]:
    """重复运行场景并返回最优耗时。每次使用独立 run_dir，避免文件已存在。"""
    best = float("inf")
    core: Optional[CorePython] = None
    for i in range(REPEATS):
        run_dir = root / f"skip={skip_store}-{i}"
        run_dir.mkdir(parents=True)
        core = _make_core(run_dir, skip_store)
        best = min(best, _write_all(core, records))
        assert core._store is not None
        core._store.close()
    assert core is not None
    return best, core


def _report(result: Dict[str, object]) -> None:
    print("\n" + "=" * 68)
    print("  core.skip_store persistence benchmark")
    print("=" * 68)
    for k, v in result.items():
        print(f"  {k:34s}: {v}")
    print("=" * 68)


# ===========================================================================
# Benchmark
# ===========================================================================


def test_bench_skip_store_vs_persist(tmp_path):
    """落盘 vs 完全跳过持久化的端到端对比（含 protobuf 序列化）。"""
    records = _make_records(NUM_KEYS, NUM_STEPS)

    normal_time, normal_core = _run_scenario(tmp_path, skip_store=False, records=records)
    skip_time, skip_core = _run_scenario(tmp_path, skip_store=True, records=records)

    # ---- 数据完整性：normal 回读全部记录；skip 不产生任何文件 ----
    reader = DataStoreReader()
    reader.open(str(normal_core._ctx.run_file))
    read_count = sum(1 for _ in reader)
    reader.close()
    skip_file_exists = skip_core._ctx.run_file.exists()

    speedup = normal_time / skip_time if skip_time > 0 else float("inf")
    result: Dict[str, object] = {
        "num_keys": NUM_KEYS,
        "num_steps": NUM_STEPS,
        "num_records": NUM_RECORDS,
        "batch_size": BATCH_SIZE,
        "persist_total_sec": round(normal_time, 4),
        "persist_rec_per_sec": round(NUM_RECORDS / normal_time, 1),
        "persist_us_per_rec": round(normal_time / NUM_RECORDS * 1e6, 3),
        "skip_total_sec": round(skip_time, 4),
        "skip_rec_per_sec": round(NUM_RECORDS / skip_time, 1) if skip_time > 0 else "inf",
        "skip_us_per_rec": round(skip_time / NUM_RECORDS * 1e6, 3),
        "speedup_x": round(speedup, 2),
        "persist_file_bytes": normal_core._ctx.run_file.stat().st_size,
        "persist_reader_count": read_count,
        "skip_file_created": skip_file_exists,
    }
    _report(result)

    # 正确性：落盘可完整回读，skip 完全不落盘
    assert read_count == NUM_RECORDS, f"Data integrity check failed: wrote {NUM_RECORDS}, read {read_count}"
    assert skip_file_exists is False, "skip_store must not create a local run file"
    # 性能：跳过序列化 + 落盘必须不慢于正常落盘
    assert skip_time <= normal_time, f"skip_store path ({skip_time:.4f}s) slower than persist ({normal_time:.4f}s)"
