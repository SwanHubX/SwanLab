"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: 指标上传线程，基于 protobuf Record
"""

from .collector import Collector
from .sender import (
    CoreTransportConfig,
    NoopRecordTransport,
    RecordTransport,
    create_record_transport,
    trace_records,
    upload_records,
)
from .uploader import Uploader

__all__ = [
    "CoreTransportConfig",
    "Collector",
    "NoopRecordTransport",
    "RecordTransport",
    "Uploader",
    "create_record_transport",
    "trace_records",
    "upload_records",
]
