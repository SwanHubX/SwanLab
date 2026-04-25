"""
@author: cunyue
@file: __init__.py
@time: 2026/4/25 16:17
@description: 列处理模块
"""

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.core_python.api.upload import UploadColumn


@safe.block(message="Failed to encode column record")
def encode(record: ColumnRecord) -> UploadColumn:
    """
    将列记录编码为后端所需的格式（DTO）
    """
    ...


@safe.block(message="Failed to decode column record")
def decode(data: bytes):
    """
    将后端返回的数据解码为列记录
    """
    raise NotImplementedError("Waiting for resume feature")
