"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: 指标上传线程，基于 protobuf Record
"""

from .thread import ThreadPool
from .upload import (
    CoreTransportConfig,
    NoopRecordTransport,
    RecordLike,
    RecordTransport,
    create_record_transport,
    load_record,
    trace_records,
    upload_records,
)

__all__ = [
    "CoreTransportConfig",
    "NoopRecordTransport",
    "RecordLike",
    "RecordTransport",
    "create_record_transport",
    "load_record",
    "trace_records",
    "ThreadPool",
    "upload_records",
]
