"""ScoreClock：分数事件时间必须按槽位严格递增，并跨进程重启保持。"""

import json

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScoreRecord
from swanlab.sdk.internal.core_python.metrics.score import ScoreClock

SLOT = ("eval0000000000000000a", "validation", "accuracy")


def test_same_millisecond_is_strictly_increasing(monkeypatch):
    clock = ScoreClock()
    # 固定墙钟，模拟同一毫秒内的连续更新
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))

    first = clock.next_timestamp(SLOT)
    second = clock.next_timestamp(SLOT)

    assert first.ToNanoseconds() == 1_000 * 10**6
    assert second.ToNanoseconds() == 1_001 * 10**6


def test_clock_rollback_keeps_increasing(monkeypatch):
    clock = ScoreClock()
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 5_000))
    first = clock.next_timestamp(SLOT)

    # 时钟回拨到更早的时间
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))
    second = clock.next_timestamp(SLOT)

    assert second.ToNanoseconds() > first.ToNanoseconds()


def test_slots_are_independent(monkeypatch):
    clock = ScoreClock()
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))

    clock.next_timestamp(SLOT)
    other = clock.next_timestamp((SLOT[0], "train", SLOT[2]))

    # 另一个槽位不受影响，仍取当前墙钟
    assert other.ToNanoseconds() == 1_000 * 10**6


def test_watermark_survives_restart(tmp_path, monkeypatch):
    state = tmp_path / "score_clock.json"
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 2_000))
    first = ScoreClock(state).next_timestamp(SLOT)

    # 新进程 + 时钟回拨，仍必须严格大于重启前的版本
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))
    second = ScoreClock(state).next_timestamp(SLOT)

    assert second.ToNanoseconds() > first.ToNanoseconds()
    assert json.loads(state.read_text(encoding="utf-8"))


def test_missing_ledger_recovers_watermark_from_records(tmp_path, monkeypatch):
    """账本丢失时用已落盘的分数记录重建水位，而不是退回墙钟。"""
    state = tmp_path / "score_clock.json"
    timestamp = Timestamp()
    timestamp.FromNanoseconds(9_000 * 10**6)
    record = ScoreRecord(eval_id=SLOT[0], dataset=SLOT[1], key=SLOT[2], timestamp=timestamp)
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))

    clock = ScoreClock(state, recover=lambda: [record])

    assert clock.next_timestamp(SLOT).ToNanoseconds() == 9_001 * 10**6


def test_corrupted_ledger_recovers_from_records(tmp_path, monkeypatch):
    state = tmp_path / "score_clock.json"
    state.write_text("not json", encoding="utf-8")
    timestamp = Timestamp()
    timestamp.FromNanoseconds(8_000 * 10**6)
    record = ScoreRecord(eval_id=SLOT[0], dataset=SLOT[1], key=SLOT[2], timestamp=timestamp)
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))

    assert ScoreClock(state, recover=lambda: [record]).next_timestamp(SLOT).ToNanoseconds() == 8_001 * 10**6


def test_stale_ledger_is_merged_with_persisted_records(tmp_path, monkeypatch):
    """账本可解析但落后（水位写盘失败）时，必须以记录里的真实最大水位为准。

    回归：1000 已落账本 → 1001 写盘失败但记录已落盘 → 重启仍读到 1000，
    只看账本会再次分配 1001，让两个不同分值撞在同一版本上。
    """
    state = tmp_path / "score_clock.json"
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))
    first = ScoreClock(state, recover=list)
    assert first.next_timestamp(SLOT).ToNanoseconds() == 1_000 * 10**6

    timestamp = Timestamp()
    timestamp.FromNanoseconds(1_001 * 10**6)
    record = ScoreRecord(eval_id=SLOT[0], dataset=SLOT[1], key=SLOT[2], timestamp=timestamp)

    # 账本仍停在 1000（那次写盘失败），但记录里已经有 1001
    restarted = ScoreClock(state, recover=lambda: [record])
    assert restarted.next_timestamp(SLOT).ToNanoseconds() == 1_002 * 10**6


def test_recovery_failure_blocks_writes_instead_of_guessing(tmp_path, monkeypatch):
    """无法确定安全水位时拒写该进程的分数，而不是猜一个可能重复的版本。"""
    state = tmp_path / "score_clock.json"
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 3_000))

    def _boom():
        raise OSError("store unreadable")

    clock = ScoreClock(state, recover=_boom)

    with pytest.raises(RuntimeError, match="refusing to assign a score version"):
        clock.next_timestamp(SLOT)


def test_unwritable_ledger_keeps_allocation_monotonic(tmp_path, monkeypatch):
    """账本写不进去时必须留痕，且同进程内仍严格递增。"""
    state = tmp_path / "missing-dir" / "score_clock.json"
    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))
    clock = ScoreClock(state)

    first = clock.next_timestamp(SLOT)
    second = clock.next_timestamp(SLOT)

    assert second.ToNanoseconds() == first.ToNanoseconds() + 10**6
    assert clock._degraded is True


def test_restore_seeds_watermark_from_records(monkeypatch):
    clock = ScoreClock()
    timestamp = Timestamp()
    timestamp.FromNanoseconds(9_000 * 10**6)
    clock.restore([ScoreRecord(eval_id=SLOT[0], dataset=SLOT[1], key=SLOT[2], timestamp=timestamp)])

    monkeypatch.setattr(ScoreClock, "_now_ms", staticmethod(lambda: 1_000))
    assert clock.next_timestamp(SLOT).ToNanoseconds() == 9_001 * 10**6
