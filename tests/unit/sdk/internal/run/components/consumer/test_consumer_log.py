"""BackgroundConsumer 指标日志分流测试：无 define 走 _log_plain，有 define 走 _log_define。

关键契约：
- 无 define：纯构建 data record，不产出 ColumnRecord（列由 core 收到数据后自动创建）；
- 有 define：define 过的 key 首次 log 物化一次列，后续 log 不再产出；
- key 首次 log 后到达的 define 无法认领该 key（automatic 钉住）。
"""

from unittest.mock import MagicMock

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus.events import MetricDefineEvent, MetricLogEvent
from swanlab.sdk.internal.pkg import console as pkg_console


def _log_event(data: dict, step: int = 1) -> MetricLogEvent:
    return MetricLogEvent(data=data, step=step, timestamp=Timestamp(seconds=step))


class TestPlainPath:
    """无 define 路径：纯构建，零列产出。"""

    def test_no_define_no_columns(self, make_consumer):
        consumer = make_consumer()
        consumer._handle_event(_log_event({"loss": 0.5}))

        assert len(consumer._scalar_batch) == 1
        assert consumer._scalar_batch[0].key == "loss"
        assert consumer._column_batch == []

    def test_forged_system_key_filtered_with_warning(self, make_consumer, monkeypatch):
        warning = MagicMock()
        monkeypatch.setattr(pkg_console, "warning", warning)
        consumer = make_consumer()
        consumer._handle_event(_log_event({"__swanlab__.cpu": 1, "loss": 2}))

        warning.assert_called_once()
        assert len(consumer._scalar_batch) == 1
        assert consumer._scalar_batch[0].key == "loss"
        assert consumer._column_batch == []


class TestDefinePath:
    """有 define 路径：define 过的 key 物化一次列。"""

    def test_define_materializes_column_once(self, make_consumer):
        consumer = make_consumer()
        consumer._handle_event(MetricDefineEvent(key="train/loss", section_name="Train"))

        consumer._handle_event(_log_event({"train/loss": 0.5}))
        consumer._handle_event(_log_event({"train/loss": 0.4}, step=2))

        assert len(consumer._column_batch) == 1
        assert consumer._column_batch[0].column_key == "train/loss"
        assert consumer._column_batch[0].section_name == "Train"
        assert len(consumer._scalar_batch) == 2

    def test_undefined_key_in_define_run_stays_columnless(self, make_consumer):
        """同一 run 内未 define 的 key 仍不产出列（core 自动建）。"""
        consumer = make_consumer()
        consumer._handle_event(MetricDefineEvent(key="train/loss", section_name="Train"))

        consumer._handle_event(_log_event({"train/loss": 0.5, "other": 1}))

        assert len(consumer._column_batch) == 1
        assert consumer._column_batch[0].column_key == "train/loss"
        assert len(consumer._scalar_batch) == 2

    def test_late_define_cannot_claim_logged_key(self, make_consumer):
        """key 首次 log 后到达的 define 无法认领该 key（automatic 钉住）。"""
        consumer = make_consumer()
        consumer._handle_event(_log_event({"loss": 0.5}))  # 无 define，先 log
        consumer._handle_event(MetricDefineEvent(key="loss", section_name="Train"))
        consumer._handle_event(_log_event({"loss": 0.4}, step=2))

        assert consumer._column_batch == []
        assert len(consumer._scalar_batch) == 2
