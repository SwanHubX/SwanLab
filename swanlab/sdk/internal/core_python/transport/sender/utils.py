"""
@author: caddiesnew
@file: utils.py
@time: 2026/3/30 15:33
@description: sender 通用工具：异常转换、可重试判定与 DTO 编码
"""

import hashlib
from typing import Optional, Sequence, Union

import requests

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import adapter
from swanlab.sdk.typings.core_python.api.upload import UploadColumn
from swanlab.sdk.typings.core_python.transport import FailureSource, FailureStage, UploadFailure


def compute_md5(source: Union[str, bytes], chunk_size: int = 8 * 1024 * 1024) -> str:
    """计算文件路径或字节数据的 MD5。"""
    md5 = hashlib.md5()
    if isinstance(source, bytes):
        md5.update(source)
    else:
        with open(source, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5.update(chunk)
    return md5.hexdigest()


def transport_failure(
    records: Sequence[Record],
    error: BaseException,
    *,
    source: FailureSource,
    stage: Optional[FailureStage] = None,
    message: Optional[str] = None,
) -> UploadFailure:
    """将上传异常转换为 UploadFailure。"""
    if isinstance(error, ApiError):
        return UploadFailure(
            source=source,
            stage=stage,
            message=message if message is not None else str(error.message),
            records=list(records),
            status=getattr(error.response, "status_code", None),
            code=error.code,
            exception=error,
        )
    return UploadFailure(
        source=source,
        stage=stage,
        message=message if message is not None else type(error).__name__,
        records=list(records),
        exception=error,
    )


def is_retryable_error(error: requests.RequestException) -> bool:
    """判断请求失败是否属于可恢复的网络或服务端故障。"""
    if not isinstance(error, ApiError):
        return True
    status = getattr(error.response, "status_code", None)
    return status is None or status >= 500


def encode_column(record: ColumnRecord) -> UploadColumn:
    """将列记录编码为后端所需的格式（DTO）。"""
    # TODO 在 v0.11 版本 精简一下 column 的类型定义
    column: UploadColumn = {"key": record.column_key, "type": adapter.column[record.column_type]}
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["type"] = "SYSTEM"
    if record.HasField("x_axis"):
        column["xAxis"] = record.x_axis
    if record.section_name:
        column["sectionName"] = record.section_name
    if record.hidden:
        column["hidden"] = True
    return column


__all__ = [
    "compute_md5",
    "encode_column",
    "is_retryable_error",
    "transport_failure",
]
