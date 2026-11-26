"""log_score：输入校验、记录构建与本地落盘（不触网）。"""

from pathlib import Path

import pytest

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScoreRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.cmd.init import init
from swanlab.sdk.internal.bus.events import ScoreLogEvent
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.run import Run
from swanlab.sdk.internal.run.components.consumer.builder import RecordBuilder

EVAL_ID = "eval0000000000000000a"


def _read_records(run: Run) -> list[Record]:
    reader = DataStoreReader()
    reader.open(str(run._ctx.run_file))
    records: list[Record] = []
    for payload in reader:
        record = Record()
        record.ParseFromString(payload)
        records.append(record)
    return records


@pytest.fixture
def local_run(tmp_path: Path):
    run = init(mode="local", log_dir=str(tmp_path / "swanlog"), reinit=True)
    try:
        yield run
    finally:
        run.finish()


def test_log_score_persists_millisecond_event_time(local_run: Run):
    local_run.log_score({"accuracy": 0.25}, eval_id=EVAL_ID, dataset="validation")
    local_run.finish()

    scores = [r.score for r in _read_records(local_run) if r.HasField("score")]
    assert len(scores) == 1
    score = scores[0]
    assert (score.eval_id, score.dataset, score.key) == (EVAL_ID, "validation", "accuracy")
    assert score.value.number == pytest.approx(0.25)
    # 毫秒精度：纳秒必须是整毫秒
    assert score.timestamp.nanos % 1_000_000 == 0


def test_same_key_versions_are_strictly_increasing(local_run: Run):
    local_run.log_score({"accuracy": 0.25}, eval_id=EVAL_ID)
    local_run.log_score({"accuracy": 0.75}, eval_id=EVAL_ID)
    local_run.finish()

    times = [r.score.timestamp.ToNanoseconds() for r in _read_records(local_run) if r.HasField("score")]
    assert len(times) == 2
    assert times[1] > times[0]


def test_invalid_eval_id_is_rejected(local_run: Run):
    local_run.log_score({"accuracy": 0.25}, eval_id="too-short")
    local_run.finish()

    assert not [r for r in _read_records(local_run) if r.HasField("score")]


def test_overall_score_dataset_is_preserved(local_run: Run):
    """空串 dataset 是"整体分数"槽位，不能被改写成 default。"""
    local_run.log_score({"accuracy": 0.5}, eval_id=EVAL_ID, dataset="")
    local_run.finish()

    scores = [r.score for r in _read_records(local_run) if r.HasField("score")]
    assert [s.dataset for s in scores] == [""]


def test_dataset_length_is_bounded(local_run: Run):
    local_run.log_score({"accuracy": 0.5}, eval_id=EVAL_ID, dataset="d" * 513)
    local_run.finish()

    assert not [r for r in _read_records(local_run) if r.HasField("score")]


def test_string_value_is_rejected_instead_of_becoming_nan(local_run: Run):
    """Scalar.transform 会把 "inf" 静默变成 NaN；分数必须在归一化之前拒绝字符串。"""
    local_run.log_score({"accuracy": "inf", "f1": 0.5}, eval_id=EVAL_ID)
    local_run.finish()

    scores = [r.score for r in _read_records(local_run) if r.HasField("score")]
    assert [s.key for s in scores] == ["f1"]


def test_infinite_score_is_skipped_not_silently_degraded(local_run: Run):
    local_run.log_score({"accuracy": float("inf"), "f1": 0.5}, eval_id=EVAL_ID)
    local_run.finish()

    scores = [r.score for r in _read_records(local_run) if r.HasField("score")]
    assert [s.key for s in scores] == ["f1"]


def test_builder_assigns_independent_timestamps_per_key(tmp_path: Path):
    run = init(mode="local", log_dir=str(tmp_path / "swanlog"), reinit=True)
    try:
        builder = RecordBuilder(run._ctx)
        records = builder.build_score(ScoreLogEvent(eval_id=EVAL_ID, dataset="validation", scores={"a": 1, "b": 2}))
        assert [r.key for r in records] == ["a", "b"]
        assert isinstance(records[0], ScoreRecord)
        # 不同槽位各自取墙钟，互不挤占
        assert records[0].timestamp.nanos % 1_000_000 == 0
    finally:
        run.finish()
