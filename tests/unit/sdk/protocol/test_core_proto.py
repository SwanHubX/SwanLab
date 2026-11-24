"""
@author: likang
@file: test_core_proto.py
@time: 2026/11/3
@description: Core proto 启动错误分类契约
"""

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
