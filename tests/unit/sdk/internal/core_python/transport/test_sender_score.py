"""分数上传：按 eval 分组、时间戳原样透传、非法记录不静默补时间。"""

from pathlib import Path
from unittest.mock import patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarValue, ScoreRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender

EVAL_A = "eval0000000000000000a"


class _FakeResponse:
    def __init__(self, status_code: int) -> None:
        self.status_code = status_code
        self.request = None
        self.url = "https://api.test/house/metrics"


EVAL_B = "eval0000000000000000b"


def _make_sender(tmp_path: Path) -> HttpRecordSender:
    ctx = CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=tmp_path,
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
        )
    )
    ctx.set_online_params(
        username="alice",
        project="demo",
        project_id="project-id",
        project_version=1,
        experiment_id="experiment-id",
    )
    return HttpRecordSender(ctx=ctx)


def _make_score(
    *,
    eval_id: str = EVAL_A,
    dataset: str = "validation",
    key: str = "accuracy",
    number: float = 0.25,
    millis: int = 1_766_000_000_001,
    with_timestamp: bool = True,
) -> Record:
    record = ScoreRecord(eval_id=eval_id, dataset=dataset, key=key, value=ScalarValue(number=number))
    if with_timestamp:
        ts = Timestamp()
        ts.FromNanoseconds(millis * 10**6)
        record.timestamp.CopyFrom(ts)
    return Record(score=record)


def test_score_batches_are_split_per_eval_id(tmp_path: Path):
    sender = _make_sender(tmp_path)
    records = [_make_score(), _make_score(eval_id=EVAL_B, key="f1", number=0.5)]

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score") as mock_upload:
        result = sender.upload_score(records)

    calls = {call.args[1]: call.kwargs["metrics"] for call in mock_upload.call_args_list}
    assert set(calls) == {EVAL_A, EVAL_B}
    assert calls[EVAL_A][0]["key"] == "accuracy"
    assert calls[EVAL_B][0]["key"] == "f1"
    assert len(result.confirmed_records) == 2
    assert result.remaining_records == []


def test_dataset_empty_string_is_not_rewritten(tmp_path: Path):
    """空串是"整体分数"槽位，不能改写成 default。"""
    sender = _make_sender(tmp_path)

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score") as mock_upload:
        sender.upload_score([_make_score(dataset="")])

    assert mock_upload.call_args.kwargs["metrics"][0]["dataset"] == ""


def test_rejected_eval_does_not_block_the_next_eval(tmp_path: Path):
    """A 返回 400 只终止 A，B 仍要发出去（评审 P1）。"""
    sender = _make_sender(tmp_path)
    records = [_make_score(), _make_score(eval_id=EVAL_B, key="f1", number=0.5)]

    def _reject_a(project_id, eval_id, *, metrics):
        if eval_id == EVAL_A:
            raise ApiError(_FakeResponse(400), method="POST", code=400, message="metrics[0].create_time")

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score", side_effect=_reject_a) as mock_upload:
        result = sender.upload_score(records)

    # 两个 eval 都被尝试过，A 的记录被确认为永久失败，B 被确认成功
    assert [call.args[1] for call in mock_upload.call_args_list] == [EVAL_A, EVAL_B]
    assert len(result.confirmed_records) == 2
    assert result.remaining_records == []
    assert [f.records[0].score.eval_id for f in result.failures] == [EVAL_A]


def test_retryable_failure_keeps_every_unsent_eval(tmp_path: Path):
    sender = _make_sender(tmp_path)
    records = [_make_score(), _make_score(eval_id=EVAL_B, key="f1", number=0.5)]

    def _server_error(project_id, eval_id, *, metrics):
        raise ApiError(_FakeResponse(503), method="POST", code=503, message="unavailable")

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score", side_effect=_server_error):
        result = sender.upload_score(records)

    # 5xx 可重试：本组和后续组都要留在 remaining，不能被丢弃
    assert result.confirmed_records == []
    assert len(result.remaining_records) == 2


def test_score_payload_keeps_millisecond_event_time(tmp_path: Path):
    sender = _make_sender(tmp_path)

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score") as mock_upload:
        sender.upload_score([_make_score()])

    metric = mock_upload.call_args.kwargs["metrics"][0]
    assert metric["dataset"] == "validation"
    assert metric["create_time"] == "2025-12-17T19:33:20.001Z"
    assert metric["data"] == 0.25


def test_nan_is_encoded_as_lowercase_nan(tmp_path: Path):
    sender = _make_sender(tmp_path)

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score") as mock_upload:
        sender.upload_score([_make_score(number=float("nan"))])

    assert mock_upload.call_args.kwargs["metrics"][0]["data"] == "nan"


def test_record_without_event_time_is_skipped(tmp_path: Path):
    """缺时间戳时不能补 now：House 会以最大时间戳覆盖合法旧版本。"""
    sender = _make_sender(tmp_path)

    record = _make_score(with_timestamp=False)
    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_score") as mock_upload:
        result = sender.upload_score([record])

    # 没有合法记录时不发请求，且该记录被标记为确定性失败（不再重试）
    assert not mock_upload.called
    assert result.confirmed_records == [record]
    assert result.failures and "event time" in result.failures[0].message
