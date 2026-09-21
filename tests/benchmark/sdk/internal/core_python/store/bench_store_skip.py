"""
@author: nexisato
@file: bench_store_skip.py
@time: 2026/9/11
@description: core.skip_store 本地持久化开销基准测试

对比 online mode 落盘（skip_store=False）与完全跳过持久化（skip_store=True）
在相同写入负载下的性能，量化“跳过 protobuf 序列化 + 本地落盘”带来的收益。

负载规模与 bench_metrics_steps.py 对齐：200 key × 5000 step
（总记录数亦与 bench_store_fsync.py 的 NUM_RECORDS 一致）；额外叠加：
  - 1,000 条 media，每条 32 KiB）
  - 100 个 save file（源文件各 32 KiB）

三类记录在两种模式下均按生产路径构造，差异即 skip_store 真正省下的本地工作：
  标量：
    - skip_store=False: 逐条 ``Record.SerializeToString()`` 后写入 LevelDB log
    - skip_store=True : 仅累加未持久化计数，不序列化
  媒体（对应 Image.transform）：
    - skip_store=False: 32 KiB 写入 ``media/image/``（``fs.safe_write`` 含 fsync），
      MediaItem 仅带 filename/sha256/size
    - skip_store=True : 字节流内联进 ``MediaItem.payload``，不落盘
  文件保存（对应 Core._handle_custom_save）：
    - skip_store=False: 在 ``files/`` 建立软链接镜像并填充 target_path
    - skip_store=True : 不建镜像，直接引用 source_path

关注指标：
  1. 各类型 persist / skip 总耗时与加速比
  2. 标量吞吐量 (rec/s)，媒体写入带宽 (MiB/s)
  3. 数据完整性（persist 标量可完整回读、媒体文件数与 save 镜像数正确；skip 不落盘）

用法：
  uv run pytest tests/benchmark/sdk/internal/core_python/store/bench_store_skip.py -v -s
"""

import hashlib
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import (
    MediaItem,
    MediaRecord,
    MediaValue,
    ScalarRecord,
    ScalarValue,
)
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord, SaveType
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.core import CorePython
from swanlab.sdk.internal.core_python.store import DataStoreReader, DataStoreWriter
from swanlab.sdk.internal.core_python.watcher import create_save_links
from swanlab.sdk.internal.pkg import adapter, fs

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

MEDIA_COUNT = 1_000
MEDIA_SIZE = 32 * 1024  # 每条媒体 32 KiB
SAVE_COUNT = 100
SAVE_SIZE = 32 * 1024  # 每个 CUSTOM save 源文件 32 KiB


# ===========================================================================
# 记录构造
# ===========================================================================


def _make_scalar_records() -> List[Record]:
    """预生成 step × key 的标量 Record，避免把构建开销算进持久化基准。

    key 命名与 bench_metrics_steps.py 的 ``k{i}`` 一致，便于横向对比。
    """
    records: List[Record] = []
    num = 0
    for step in range(1, NUM_STEPS + 1):
        for k in range(NUM_KEYS):
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


def _make_media_records(skip_store: bool) -> Tuple[List[Record], List[Tuple[str, bytes]]]:
    """构造 1,000 条 32 KiB 的 Image media。

    返回 (records, files)；``files`` 仅在 skip_store=False 时非空，为需落盘的
    (filename, content)。skip_store=True 时 content 内联进 MediaItem.payload。
    """
    records: List[Record] = []
    files: List[Tuple[str, bytes]] = []
    for i in range(MEDIA_COUNT):
        content = bytes([i % 256]) * MEDIA_SIZE
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{i:04d}-{sha256[:8]}.png"
        item = MediaItem(filename=filename, sha256=sha256, size=len(content), caption="")
        if skip_store:
            item.payload = content
        else:
            files.append((filename, content))
        records.append(
            Record(
                num=NUM_RECORDS + i + 1,
                media=MediaRecord(
                    key=f"media/img_{i % 10}",
                    step=i + 1,
                    type=ColumnType.COLUMN_TYPE_IMAGE,
                    value=MediaValue(items=[item]),
                ),
            )
        )
    return records, files


def _make_save_entries(source_paths: List[Path]) -> List[SaveRecord]:
    """构造 CUSTOM SaveRecord；CUSTOM 的 payload 恒空，只引用 source_path。"""
    return [
        SaveRecord(
            name=f"checkpoints/model_{i:03d}.bin",
            source_path=str(src),
            policy=SavePolicy.SAVE_POLICY_END,
            type=SaveType.SAVE_TYPE_CUSTOM,
        )
        for i, src in enumerate(source_paths)
    ]


def _to_save_envelopes(saves: List[SaveRecord]) -> List[Record]:
    """将 SaveRecord 包成 Record envelope（target_path 由 create_save_links 填充后再拷贝）。"""
    total = NUM_RECORDS + MEDIA_COUNT
    return [Record(num=total + i + 1, save=save) for i, save in enumerate(saves)]


def _make_source_files(root: Path) -> List[Path]:
    """创建 SAVE_COUNT 个 32 KiB 源文件（基准计时之外，仅作为 save 的引用目标）。"""
    src_dir = root / "sources"
    src_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for i in range(SAVE_COUNT):
        path = src_dir / f"model_{i:03d}.bin"
        path.write_bytes(bytes([i % 256]) * SAVE_SIZE)
        paths.append(path)
    return paths


# ===========================================================================
# 持久化执行
# ===========================================================================


def _make_core(root: Path, skip_store: bool, tag: str) -> CorePython:
    """构造仅用于持久化路径的 CorePython（不启动 transport / heartbeat）。"""
    run_dir = root / f"{tag}-skip={skip_store}"
    run_dir.mkdir(parents=True)
    core = CorePython(mode="online")
    core._ctx = CoreContext(
        config=CoreConfig(
            run_id=f"bench-{tag}",
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


def _write_records(core: CorePython, records: List[Record]) -> None:
    """按 BATCH_SIZE 分批调用生产路径 _store_records。"""
    for offset in range(0, len(records), BATCH_SIZE):
        core._store_records(records[offset : offset + BATCH_SIZE])


def _run_scalar(root: Path, skip_store: bool) -> Tuple[float, CorePython]:
    records = _make_scalar_records()
    best = float("inf")
    core: Optional[CorePython] = None
    for i in range(REPEATS):
        core = _make_core(root, skip_store, f"scalar-{i}")
        start = time.perf_counter()
        _write_records(core, records)
        best = min(best, time.perf_counter() - start)
        assert core._store is not None
        core._store.close()
    assert core is not None
    return best, core


def _run_media(root: Path, skip_store: bool) -> Tuple[float, CorePython]:
    best = float("inf")
    core: Optional[CorePython] = None
    for i in range(REPEATS):
        core = _make_core(root, skip_store, f"media-{i}")
        records, files = _make_media_records(skip_store)
        start = time.perf_counter()
        if not skip_store:
            media_dir = core._ctx.media_dir / adapter.medium[ColumnType.COLUMN_TYPE_IMAGE]
            fs.safe_mkdir(media_dir)
            for filename, content in files:
                fs.safe_write(media_dir / filename, content, mode="wb")
        _write_records(core, records)
        best = min(best, time.perf_counter() - start)
        assert core._store is not None
        core._store.close()
    assert core is not None
    return best, core


def _run_saves(root: Path, skip_store: bool, source_paths: List[Path]) -> Tuple[float, CorePython]:
    best = float("inf")
    core: Optional[CorePython] = None
    for i in range(REPEATS):
        core = _make_core(root, skip_store, f"save-{i}")
        saves = _make_save_entries(source_paths)
        start = time.perf_counter()
        if not skip_store:
            create_save_links(saves, core._ctx.files_dir)
        _write_records(core, _to_save_envelopes(saves))
        best = min(best, time.perf_counter() - start)
        assert core._store is not None
        core._store.close()
    assert core is not None
    return best, core


# ===========================================================================
# 报告
# ===========================================================================


def _report(rows: List[Tuple[str, str, float, float]], notes: Dict[str, object]) -> None:
    print("\n" + "=" * 78)
    print("  core.skip_store persistence benchmark")
    print("=" * 78)
    for k, v in notes.items():
        print(f"  {k:16s}: {v}")
    print("-" * 78)
    print(f"  {'category':<10}{'scale':<16}{'persist(s)':>14}{'skip(s)':>14}{'speedup':>12}")
    print("-" * 78)
    for name, scale, persist, skip in rows:
        speedup = f"{persist / skip:,.1f}x" if skip > 0 else "inf"
        print(f"  {name:<10}{scale:<16}{persist:>14.4f}{skip:>14.6f}{speedup:>12}")
    total_persist = sum(r[2] for r in rows)
    total_skip = sum(r[3] for r in rows)
    total_speedup = f"{total_persist / total_skip:,.1f}x" if total_skip > 0 else "inf"
    print("-" * 78)
    print(f"  {'TOTAL':<10}{'':<16}{total_persist:>14.4f}{total_skip:>14.6f}{total_speedup:>12}")
    print("=" * 78)


# ===========================================================================
# Benchmark
# ===========================================================================


def test_bench_skip_store_vs_persist(tmp_path):
    """标量 + 媒体 + 文件保存：落盘 vs 完全跳过持久化的端到端对比。"""
    source_paths = _make_source_files(tmp_path)

    p_scalar, p_scalar_core = _run_scalar(tmp_path, skip_store=False)
    s_scalar, s_scalar_core = _run_scalar(tmp_path, skip_store=True)
    p_media, p_media_core = _run_media(tmp_path, skip_store=False)
    s_media, s_media_core = _run_media(tmp_path, skip_store=True)
    p_save, p_save_core = _run_saves(tmp_path, skip_store=False, source_paths=source_paths)
    s_save, s_save_core = _run_saves(tmp_path, skip_store=True, source_paths=source_paths)

    # ---- 数据完整性：persist 全部落盘，skip 全部不落盘 ----
    reader = DataStoreReader()
    reader.open(str(p_scalar_core._ctx.run_file))
    scalar_read_count = sum(1 for _ in reader)
    reader.close()

    media_dir = p_media_core._ctx.media_dir / adapter.medium[ColumnType.COLUMN_TYPE_IMAGE]
    media_files = list(media_dir.glob("*.png"))
    media_bytes = sum(f.stat().st_size for f in media_files)

    save_links = [p for p in p_save_core._ctx.files_dir.rglob("*") if not p.is_dir()]
    skip_files_absent = all(not core._ctx.run_file.exists() for core in (s_scalar_core, s_media_core, s_save_core))

    rows: List[Tuple[str, str, float, float]] = [
        ("scalars", f"{NUM_KEYS}x{NUM_STEPS}", p_scalar, s_scalar),
        ("media", f"{MEDIA_COUNT}x{MEDIA_SIZE}B", p_media, s_media),
        ("saves", f"{SAVE_COUNT} files", p_save, s_save),
    ]
    notes: Dict[str, object] = {
        "records": f"{NUM_RECORDS} scalars + {MEDIA_COUNT} media + {SAVE_COUNT} saves",
        "persist_scalar_throughput": f"{NUM_RECORDS / p_scalar:,.0f} rec/s",
        "persist_media_bandwidth": f"{media_bytes / p_media / 1024 / 1024:,.1f} MiB/s",
        "persist_run_file": f"{p_scalar_core._ctx.run_file.stat().st_size:,} bytes",
        "persist_media_files": f"{len(media_files)} files / {media_bytes:,} bytes",
        "persist_save_links": f"{len(save_links)} links",
        "skip_local_files_absent": skip_files_absent,
    }
    _report(rows, notes)

    # 正确性断言
    assert scalar_read_count == NUM_RECORDS, f"scalar integrity failed: {scalar_read_count} != {NUM_RECORDS}"
    assert len(media_files) == MEDIA_COUNT, f"media integrity failed: {len(media_files)} files"
    assert media_bytes == MEDIA_COUNT * MEDIA_SIZE, f"media bytes mismatch: {media_bytes}"
    assert len(save_links) == SAVE_COUNT, f"save links mismatch: {len(save_links)}"
    assert skip_files_absent, "skip_store must not create a local run file"
    for name, _, persist, skip in rows:
        assert skip <= persist, f"skip_store path slower than persist for {name}: {skip:.4f}s > {persist:.4f}s"
