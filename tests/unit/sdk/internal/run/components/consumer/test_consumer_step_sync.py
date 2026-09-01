"""BackgroundConsumer step_sync 注入路径测试。"""

import math
import queue
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.sdk.internal.bus.events import MetricDefineEvent, MetricLogEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer
from swanlab.sdk.internal.run.transforms import Text


def _make_consumer(tmp_path: Path) -> BackgroundConsumer:
    ctx = SimpleNamespace(
        core=MagicMock(),
        callbacker=MagicMock(),
        media_dir=tmp_path,
        config=SimpleNamespace(settings=SimpleNamespace(core=SimpleNamespace(section_rule=0))),
    )
    run_ctx = cast(RunContext, cast(object, ctx))
    return BackgroundConsumer(run_ctx, queue.Queue())


def _timestamp(seconds: int) -> Timestamp:
    return Timestamp(seconds=seconds)


def _log(consumer: BackgroundConsumer, data: dict, step: int, seconds: int) -> None:
    consumer._handle_event(MetricLogEvent(data=data, step=step, timestamp=_timestamp(seconds)))


def _records_by_key(consumer: BackgroundConsumer, start: int = 0):
    records = consumer._scalar_batch[start:]
    return {key: [record for record in records if record.key == key] for key in {record.key for record in records}}


def test_step_sync_injects_cached_x_with_current_y_step_and_timestamp(tmp_path: Path):
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/loss": 0.8}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"epoch", "train/loss"}
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].value.number == 3
    assert records["epoch"][0].step == 1
    assert records["epoch"][0].timestamp == _timestamp(20)
    assert records["train/loss"][0].value.number == pytest.approx(0.8)
    assert records["train/loss"][0].step == 1
    assert records["train/loss"][0].timestamp == _timestamp(20)
    # 只有 define 过的 train/loss 产出列；epoch 未 define，列由 core 收到数据后自动创建
    assert [column.column_key for column in consumer._column_batch] == ["train/loss"]


@pytest.mark.parametrize(
    "data",
    [
        {"epoch": 3, "train/loss": 0.8},
        {"train/loss": 0.8, "epoch": 3},
    ],
    ids=["x-before-y", "y-before-x"],
)
def test_step_sync_does_not_inject_when_event_contains_x_regardless_of_order(tmp_path: Path, data: dict):
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))

    _log(consumer, data, step=1, seconds=20)

    records = _records_by_key(consumer)
    assert set(records) == {"epoch", "train/loss"}
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].value.number == 3


def test_step_sync_injects_shared_x_only_once_for_multiple_y_in_same_event(tmp_path: Path):
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch"))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/acc": 0.9, "train/loss": 0.1}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"epoch", "train/acc", "train/loss"}
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].value.number == 3


def test_step_sync_does_not_inject_shared_x_twice_at_same_step(tmp_path: Path):
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch"))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/acc": 0.9}, step=1, seconds=20)
    _log(consumer, {"train/loss": 0.1}, step=1, seconds=21)

    records = _records_by_key(consumer, start)
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].timestamp == _timestamp(20)


def test_step_sync_injects_default_zero_before_first_real_x(tmp_path: Path):
    """X 从未 log 过时，绑定 Y 触发注入默认 0，X 序列从首个绑定 Y 的 step 出现（对齐 wandb）。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch"))
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/acc": 0.9}, step=0, seconds=10)

    records = _records_by_key(consumer, start)
    assert set(records) == {"train/acc", "epoch"}
    assert records["epoch"][0].value.number == 0
    assert records["epoch"][0].step == 0
    assert records["epoch"][0].timestamp == _timestamp(10)


def test_step_sync_default_zero_at_each_y_step_then_forward_fill_after_real_x(tmp_path: Path):
    """无值窗口内每个 Y event 在各自 step 注入 0；真实 X 到达后恢复取最近真实值。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch"))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))

    _log(consumer, {"train/acc": 0.9}, step=0, seconds=10)
    _log(consumer, {"train/loss": 0.1}, step=1, seconds=20)
    _log(consumer, {"epoch": 3}, step=2, seconds=30)
    _log(consumer, {"train/acc": 0.8}, step=3, seconds=40)

    epoch_records = [r for r in consumer._scalar_batch if r.key == "epoch"]
    assert [(r.step, r.value.number) for r in epoch_records] == [(0, 0), (1, 0), (2, 3), (3, 3)]


def test_real_x_after_injected_x_warns_once_and_updates_next_step_cache(tmp_path: Path, monkeypatch):
    warnings = []
    monkeypatch.setattr(
        "swanlab.sdk.internal.run.components.consumer.resolver.console.warning",
        lambda message: warnings.append(message),
    )
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    _log(consumer, {"train/loss": 0.8}, step=1, seconds=20)

    _log(consumer, {"epoch": 4}, step=1, seconds=21)
    _log(consumer, {"epoch": 4}, step=1, seconds=22)
    start = len(consumer._scalar_batch)
    _log(consumer, {"train/loss": 0.7}, step=2, seconds=30)

    records = _records_by_key(consumer, start)
    assert len(warnings) == 1
    assert "real value arrives after injected value" in warnings[0]
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].value.number == 4
    assert records["epoch"][0].step == 2
    assert records["epoch"][0].timestamp == _timestamp(30)


# ============================================================
# custom X first-writer-wins：同一 X 值上首次 Y 值为准
# ============================================================


def test_duplicate_x_drops_second_y(tmp_path: Path):
    """X 只 log 一次，两次 Y log → 第二次 Y 丢弃（注入相同 X 值）。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1"))
    _log(consumer, {"x1": 5}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"acc": 0.9}, step=1, seconds=20)  # inject x1=5.0 → accepted
    _log(consumer, {"acc": 1.0}, step=2, seconds=30)  # inject x1=5.0 → rejected

    acc_records = [r for r in consumer._scalar_batch[start:] if r.key == "acc"]
    assert len(acc_records) == 1
    assert acc_records[0].value.number == pytest.approx(0.9)


def test_different_x_values_accepted(tmp_path: Path):
    """X 更新后，新 X 值的 Y 被接受。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1"))
    _log(consumer, {"x1": 5}, step=0, seconds=10)
    _log(consumer, {"acc": 0.9}, step=1, seconds=20)  # inject x1=5.0 → accepted
    _log(consumer, {"x1": 6}, step=2, seconds=30)
    start = len(consumer._scalar_batch)

    _log(consumer, {"acc": 1.0}, step=3, seconds=40)  # inject x1=6.0 → accepted

    acc_records = [r for r in consumer._scalar_batch[start:] if r.key == "acc"]
    assert len(acc_records) == 1
    assert acc_records[0].value.number == pytest.approx(1.0)


def test_explicit_same_x_value_rejected(tmp_path: Path):
    """显式 log X=5 在不同 step，Y 仍按 X 值去重。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1"))

    _log(consumer, {"x1": 5, "acc": 0.9}, step=0, seconds=10)  # accepted
    _log(consumer, {"x1": 5, "acc": 1.0}, step=1, seconds=20)  # rejected (same X value)

    acc_records = [r for r in consumer._scalar_batch if r.key == "acc"]
    assert len(acc_records) == 1
    assert acc_records[0].value.number == pytest.approx(0.9)


def test_multiple_y_sharing_x_independent(tmp_path: Path):
    """loss 和 acc 共享 x1，各自独立去重。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("loss", x_axis="x1"))
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1"))
    _log(consumer, {"x1": 5}, step=0, seconds=10)
    _log(consumer, {"loss": 0.5, "acc": 0.9}, step=1, seconds=20)  # both accepted
    start = len(consumer._scalar_batch)

    _log(consumer, {"loss": 0.4, "acc": 0.8}, step=2, seconds=30)  # both rejected

    new_records = consumer._scalar_batch[start:]
    loss_records = [r for r in new_records if r.key == "loss"]
    acc_records = [r for r in new_records if r.key == "acc"]
    assert len(loss_records) == 0
    assert len(acc_records) == 0


def test_no_custom_x_no_dedup(tmp_path: Path):
    """无 custom X 轴时，不触发去重。"""
    consumer = _make_consumer(tmp_path)

    _log(consumer, {"acc": 0.9}, step=0, seconds=10)
    _log(consumer, {"acc": 1.0}, step=1, seconds=20)

    acc_records = [r for r in consumer._scalar_batch if r.key == "acc"]
    assert len(acc_records) == 2


# ============================================================
# M1：非有限 X 值不应被注入的旧值覆盖
# 注入跳过条件必须看 event.data（用户是否提供过该 key），而非 explicit_scalars
# （后者已被 isfinite 过滤）；否则 {**data, **injected} 会用旧值覆盖用户的 nan
# ============================================================


def test_nan_x_value_not_overwritten_by_injection(tmp_path: Path):
    """X 为 NaN 时，落盘的应是 nan（用户显式值），而非被注入的旧值 3.0。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    _log(consumer, {"epoch": float("nan"), "loss": 0.8}, step=1, seconds=20)

    epoch_records = [r for r in consumer._scalar_batch if r.key == "epoch"]
    # 第二条 epoch 必须是用户显式的 nan，而非注入的 3.0
    assert math.isnan(epoch_records[-1].value.number)
    # 不应出现 step=1 处 value=3 的注入孤儿点
    injected_orphan = [r for r in epoch_records if r.step == 1 and r.value.number == 3.0]
    assert injected_orphan == []


# ============================================================
# S2：Y 被 X 去重丢弃时不留孤儿注入点、不误标 INJECTED
# 注入 commit（try_inject_x + record）必须推迟到确认至少一个存活 Y 之后
# ============================================================


def test_dropped_y_leaves_no_orphan_injected_x(tmp_path: Path):
    """Y 被 X 值去重丢弃时，不为其注入孤儿 X 点（注入只对存活 Y commit）。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("loss", x_axis="epoch"))
    _log(consumer, {"epoch": 5}, step=0, seconds=10)
    _log(consumer, {"loss": 0.9}, step=1, seconds=20)  # inject epoch=5@1，Y 存活
    start = len(consumer._scalar_batch)

    _log(consumer, {"loss": 0.8}, step=2, seconds=30)  # epoch=5 连续重复 → Y 丢弃

    new_records = consumer._scalar_batch[start:]
    # loss 被去重丢弃
    assert [r for r in new_records if r.key == "loss"] == []
    # 未为被丢的 Y 注入孤儿 epoch 点（step=2 处无 epoch）
    assert [r for r in new_records if r.key == "epoch" and r.step == 2] == []


def test_dropped_y_no_spurious_real_conflict_warning(tmp_path: Path, monkeypatch):
    """Y 被丢弃后未标 INJECTED，真实 X 晚到同 step 不触发误 REAL_CONFLICT。"""
    warnings = []
    monkeypatch.setattr(
        "swanlab.sdk.internal.run.components.consumer.resolver.console.warning",
        lambda message: warnings.append(message),
    )
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("loss", x_axis="epoch"))
    _log(consumer, {"epoch": 5}, step=0, seconds=10)
    _log(consumer, {"loss": 0.9}, step=1, seconds=20)  # inject epoch=5@1（Y 存活，标 INJECTED@1）
    _log(consumer, {"loss": 0.8}, step=2, seconds=30)  # epoch=5 重复 → Y 丢弃，不标 INJECTED@2
    _log(consumer, {"epoch": 6}, step=2, seconds=31)  # 真实 X@2，未标 INJECTED@2 → 无误告警

    conflict = [w for w in warnings if "real value arrives after injected value" in w]
    assert conflict == []


def test_step_sync_false_does_not_inject_cached_x(tmp_path: Path):
    """step_sync=False：X/Y 分次 log 时不注入 X。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch", step_sync=False))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/loss": 0.8}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"train/loss"}
    assert records["train/loss"][0].value.number == pytest.approx(0.8)
    assert records["train/loss"][0].step == 1


def test_step_sync_false_does_not_drop_y_from_cached_x(tmp_path: Path):
    """False 不用跨 step cache 去重：event 不含 X 时 Y 保留。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1", step_sync=False))
    _log(consumer, {"x1": 5, "acc": 0.9}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"acc": 1.0}, step=1, seconds=20)

    acc_records = [r for r in consumer._scalar_batch[start:] if r.key == "acc"]
    assert len(acc_records) == 1
    assert acc_records[0].value.number == pytest.approx(1.0)
    assert [r for r in consumer._scalar_batch[start:] if r.key == "x1"] == []


def test_step_sync_false_still_drops_y_when_event_has_duplicate_x(tmp_path: Path):
    """False 仍按本 event 显式 X 去重。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1", step_sync=False))
    _log(consumer, {"x1": 5, "acc": 0.9}, step=0, seconds=10)
    _log(consumer, {"x1": 5, "acc": 1.0}, step=1, seconds=20)

    acc_records = [r for r in consumer._scalar_batch if r.key == "acc"]
    assert len(acc_records) == 1
    assert acc_records[0].value.number == pytest.approx(0.9)


def test_step_sync_false_same_step_real_x_does_not_inject_or_drop(tmp_path: Path):
    """False 不用本 step 真实 X：分次 log 时 Y 保留、不注入。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("acc", x_axis="x1", step_sync=False))
    _log(consumer, {"x1": 5}, step=1, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"acc": 0.9}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"acc"}
    assert records["acc"][0].value.number == pytest.approx(0.9)


def test_step_sync_false_sibling_does_not_block_true_sibling_inject(tmp_path: Path):
    """False 的 Y 不触发注入；同 event 的 True 兄弟仍注入共享 X。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch", step_sync=False))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/acc": 0.9, "train/loss": 0.1}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"epoch", "train/acc", "train/loss"}
    assert len(records["epoch"]) == 1
    assert records["epoch"][0].value.number == 3
    assert records["epoch"][0].step == 1


def test_step_sync_false_alone_does_not_inject_even_with_true_sibling_defined(tmp_path: Path):
    """仅 False 的 Y 出现时不注入，即使存在 step_sync=True 的兄弟 rule。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/acc", x_axis="epoch", step_sync=False))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch"))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/acc": 0.9}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert set(records) == {"train/acc"}


def test_overwrite_resets_step_sync_to_true_before_first_log(tmp_path: Path):
    """overwrite=True 且未指定 step_sync 时恢复默认注入。"""
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch", step_sync=False))
    consumer._handle_event(MetricDefineEvent("train/loss", x_axis="epoch", overwrite=True))
    _log(consumer, {"epoch": 3}, step=0, seconds=10)
    start = len(consumer._scalar_batch)

    _log(consumer, {"train/loss": 0.8}, step=1, seconds=20)

    records = _records_by_key(consumer, start)
    assert "epoch" in records
    assert records["epoch"][0].value.number == 3
    assert records["epoch"][0].step == 1


def test_self_reference_x_axis_keeps_first_log_wins(tmp_path: Path):
    """自引用 x_axis（x_axis == key）：正常落盘自身值，且 first-writer-wins 语义不变。

    - 自引用 key 恒随 Y 出现，不触发 step_sync 注入；
    - 连续重复 X 值按既有 epsilon 规则丢弃该 Y 点；
    - 二次 define 不改变已冻结 concrete（仍自引用）。
    """
    consumer = _make_consumer(tmp_path)
    consumer._handle_event(MetricDefineEvent("epoch", x_axis="epoch"))

    _log(consumer, {"epoch": 1}, step=0, seconds=10)
    _log(consumer, {"epoch": 2}, step=1, seconds=20)
    _log(consumer, {"epoch": 2}, step=2, seconds=30)  # 连续重复 X → 丢弃
    _log(consumer, {"epoch": 3}, step=3, seconds=40)

    records = _records_by_key(consumer)
    assert set(records) == {"epoch"}
    assert [r.value.number for r in records["epoch"]] == [1, 2, 3]
    assert [r.step for r in records["epoch"]] == [0, 1, 3]

    # ColumnRecord 仅一条，x_axis 为自身
    columns = [c for c in consumer._column_batch if c.column_key == "epoch"]
    assert len(columns) == 1
    assert columns[0].x_axis == "epoch"

    # 二次 define（first-writer-wins：concrete 已冻结，不再产出 ColumnRecord）
    consumer._handle_event(MetricDefineEvent("epoch", x_axis="_step"))
    concrete = consumer._resolver.resolve_concrete("epoch", "SCALAR")
    assert concrete.effective.x_axis == "epoch"
    assert [c for c in consumer._column_batch if c.column_key == "epoch"] == [columns[0]]


# ============================================================
# 快路径：无 custom X 规则时预扫描整体跳过，
# scalar/media 混合 log 仍正常走 builder 路径（media 落盘、列定义产出）
# ============================================================


def test_no_custom_x_mixed_media_and_scalar_log(tmp_path: Path):
    consumer = _make_consumer(tmp_path)
    assert consumer._resolver.has_custom_x is False

    _log(consumer, {"train/loss": 0.8, "train/text": Text("hello")}, step=0, seconds=10)

    scalars = [r for r in consumer._scalar_batch if r.key == "train/loss"]
    assert len(scalars) == 1
    assert scalars[0].value.number == pytest.approx(0.8)
    media = [r for r in consumer._media_batch if r.key == "train/text"]
    assert len(media) == 1
    assert media[0].type == ColumnType.COLUMN_TYPE_TEXT
    assert media[0].value.items[0].filename.endswith(".txt")
    # 无 define：SDK 不产出 ColumnRecord，列由 core 收到数据后自动创建
    assert consumer._column_batch == []
