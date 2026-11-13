"""
@author: likang
@file: test_core_proto.py
@time: 2026/11/3
@description: Core proto 启动错误分类契约
"""

from swanlab.proto.swanlab.grpc.core.v1 import core_pb2


def test_core_error_kind_enum_and_field():
    kind = core_pb2.CORE_ERROR_KIND_TRANSIENT

    resp = core_pb2.DeliverRunStartResponse(error_kind=kind)

    assert resp.error_kind == kind
    # 9 已被 session_id 占用（Go core 的 run 会话路由依赖），error_kind 顺延到 10
    assert core_pb2.DeliverRunStartResponse.SESSION_ID_FIELD_NUMBER == 9
    assert core_pb2.DeliverRunStartResponse.ERROR_KIND_FIELD_NUMBER == 10
    assert core_pb2.DeliverRunStartResponse(success=True).error_kind == core_pb2.CORE_ERROR_KIND_UNSPECIFIED
