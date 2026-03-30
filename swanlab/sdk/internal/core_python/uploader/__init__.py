"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: 指标上传线程，基于 protobuf Record
"""

from .collector import Collector
from .sender import (
    HttpRecordTransport,
    create_record_transport,
    trace_records,
    upload_records,
)
from .thread import ThreadPool

__all__ = [
    "Collector",
    "HttpRecordTransport",
    "ThreadPool",
    "create_record_transport",
    "trace_records",
    "upload_records",
]
