"""
@author: likang
@file: test_core_proto.py
@time: 2026/11/3
@description: Core proto 启动错误分类契约
"""

import pytest

from swanlab.proto.swanlab.grpc.core.v1 import core_pb2
from swanlab.proto.swanlab.operation.v1 import operation_pb2
from swanlab.proto.swanlab.settings.core.v1 import core_pb2 as settings_core_pb2


def test_core_error_kind_enum_and_field():
    kind = core_pb2.CORE_ERROR_KIND_TRANSIENT

    resp = core_pb2.DeliverRunStartResponse(error_kind=kind)

    assert resp.error_kind == kind
    # 9 已被 session_id 占用（Go core 的 run 会话路由依赖），error_kind 顺延到 10
    assert core_pb2.DeliverRunStartResponse.SESSION_ID_FIELD_NUMBER == 9
    assert core_pb2.DeliverRunStartResponse.ERROR_KIND_FIELD_NUMBER == 10
    assert core_pb2.DeliverRunStartResponse(success=True).error_kind == core_pb2.CORE_ERROR_KIND_UNSPECIFIED


def test_core_settings_recover_online_field():
    assert settings_core_pb2.CoreSettings.RECOVER_ONLINE_FIELD_NUMBER == 21
    assert settings_core_pb2.CoreSettings().recover_online is False

    encoded = settings_core_pb2.CoreSettings(run_id="r1", recover_online=True).SerializeToString()
    decoded = settings_core_pb2.CoreSettings()
    decoded.ParseFromString(encoded)
    assert decoded.run_id == "r1"
    assert decoded.recover_online is True

    # 向后兼容：缺少字段 21 的旧消息按默认值 False 解析
    legacy_encoded = settings_core_pb2.CoreSettings(run_id="legacy").SerializeToString()
    legacy_decoded = settings_core_pb2.CoreSettings()
    legacy_decoded.ParseFromString(legacy_encoded)
    assert legacy_decoded.run_id == "legacy"
    assert legacy_decoded.recover_online is False


def test_operation_stats_cloud_status_fields():
    assert operation_pb2.CLOUD_STATE_UNSPECIFIED == 0
    assert operation_pb2.CLOUD_STATE_PENDING == 1
    assert operation_pb2.CLOUD_STATE_ATTACHED == 2
    assert operation_pb2.CLOUD_STATE_FAILED == 3

    assert operation_pb2.OperationStats.CLOUD_FIELD_NUMBER == 11
    assert operation_pb2.CloudStatus.STATE_FIELD_NUMBER == 1
    assert operation_pb2.CloudStatus.PATH_FIELD_NUMBER == 2
    assert operation_pb2.CloudStatus.ERROR_FIELD_NUMBER == 3

    empty = operation_pb2.OperationStats()
    assert not empty.HasField("cloud")
    assert empty.cloud.state == operation_pb2.CLOUD_STATE_UNSPECIFIED
    assert empty.cloud.path == ""
    assert empty.cloud.error == ""

    stats = operation_pb2.OperationStats(
        state=operation_pb2.CORE_STATE_RUNNING,
        cloud=operation_pb2.CloudStatus(
            state=operation_pb2.CLOUD_STATE_ATTACHED,
            path="/u/proj/runs/exp1",
            error="",
        ),
    )
    decoded = operation_pb2.OperationStats()
    decoded.ParseFromString(stats.SerializeToString())
    assert decoded.HasField("cloud")
    assert decoded.cloud.state == operation_pb2.CLOUD_STATE_ATTACHED
    assert decoded.cloud.path == "/u/proj/runs/exp1"
    assert decoded.cloud.error == ""

    failed = operation_pb2.OperationStats(
        cloud=operation_pb2.CloudStatus(
            state=operation_pb2.CLOUD_STATE_FAILED,
            error="401 unauthorized",
        ),
    )
    decoded_failed = operation_pb2.OperationStats()
    decoded_failed.ParseFromString(failed.SerializeToString())
    assert decoded_failed.HasField("cloud")
    assert decoded_failed.cloud.state == operation_pb2.CLOUD_STATE_FAILED
    assert decoded_failed.cloud.error == "401 unauthorized"

    # 向后兼容：缺少字段 11 的旧消息按零值解析
    legacy = operation_pb2.OperationStats(total_records=5).SerializeToString()
    legacy_decoded = operation_pb2.OperationStats()
    legacy_decoded.ParseFromString(legacy)
    assert legacy_decoded.total_records == 5
    assert not legacy_decoded.HasField("cloud")
    assert legacy_decoded.cloud.state == operation_pb2.CLOUD_STATE_UNSPECIFIED
    assert legacy_decoded.cloud.path == ""
    assert legacy_decoded.cloud.error == ""


def _make_minimal_enabled_protocol(mode):
    from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
        ConfirmRunFinishResponse,
        DeliverRunFinishRequest,
        DeliverRunFinishResponse,
        DeliverRunStartRequest,
        DeliverRunStartResponse,
        GetOperationStatsResponse,
    )
    from swanlab.sdk.protocol.core import CoreProtocol

    calls = []

    class MinimalEnabledCore(CoreProtocol):
        def _start_when_enabled(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
            calls.append(("start", start_request))
            return DeliverRunStartResponse(success=True)

        def _upsert_columns_when_enabled(self, columns) -> None:
            calls.append(("columns", columns))

        def _upsert_scalars_when_enabled(self, scalars) -> None:
            calls.append(("scalars", scalars))

        def _upsert_media_when_enabled(self, media) -> None:
            calls.append(("media", media))

        def _upsert_histograms_when_enabled(self, histograms) -> None:
            calls.append(("histograms", histograms))

        def _upsert_logs_when_enabled(self, logs) -> None:
            calls.append(("logs", logs))

        def _upsert_saves_when_enabled(self, saves) -> None:
            calls.append(("saves", saves))

        def _finish_when_enabled(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
            calls.append(("finish", finish_request))
            return DeliverRunFinishResponse(success=True)

        def _get_operation_when_online(self) -> GetOperationStatsResponse:
            return GetOperationStatsResponse(success=True)

        def _confirm_finish_when_enabled(self) -> ConfirmRunFinishResponse:
            return ConfirmRunFinishResponse(success=True)

        def fork(self):
            return self

    return MinimalEnabledCore(mode), calls


@pytest.mark.parametrize("mode", ["local", "offline", "online"])
@pytest.mark.parametrize(
    "method_name, expected_tag",
    [
        ("upsert_columns", "columns"),
        ("upsert_scalars", "scalars"),
        ("upsert_media", "media"),
        ("upsert_histograms", "histograms"),
        ("upsert_logs", "logs"),
        ("upsert_saves", "saves"),
    ],
)
def test_core_protocol_default_upsert_delegates_to_when_enabled(mode, method_name, expected_tag):
    core, calls = _make_minimal_enabled_protocol(mode)
    payload = [object()]

    getattr(core, method_name)(payload)

    assert calls == [(expected_tag, payload)]


def test_core_protocol_disabled_mode_skips_when_enabled_hooks():
    from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import DeliverRunFinishRequest, DeliverRunStartRequest

    core, calls = _make_minimal_enabled_protocol("disabled")

    assert core.deliver_run_start(DeliverRunStartRequest()).success is True
    core.upsert_columns([])
    core.upsert_scalars([])
    core.upsert_media([])
    core.upsert_histograms([])
    core.upsert_logs([])
    core.upsert_saves([])
    assert core.deliver_run_finish(DeliverRunFinishRequest()).success is True
    assert calls == []


@pytest.mark.parametrize("mode", ["local", "offline", "online"])
def test_core_protocol_default_lifecycle_delegates_to_when_enabled(mode):
    from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import DeliverRunFinishRequest, DeliverRunStartRequest

    core, calls = _make_minimal_enabled_protocol(mode)
    start_req = DeliverRunStartRequest()
    finish_req = DeliverRunFinishRequest()

    assert core.deliver_run_start(start_req).success is True
    assert core.deliver_run_finish(finish_req).success is True
    assert calls == [("start", start_req), ("finish", finish_req)]
