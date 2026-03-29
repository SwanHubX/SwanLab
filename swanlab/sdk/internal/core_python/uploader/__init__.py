"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: 基于 protobuf Record 的 HTTP batch uploader
"""

from .helper import RecordQueue
from .sender import (
    HttpRecordSender,
    NoopHttpRecordSender,
    create_http_record_sender,
    group_records_by_type,
)
from .uploader import HttpBatchUploader

__all__ = [
    "HttpBatchUploader",
    "HttpRecordSender",
    "NoopHttpRecordSender",
    "RecordQueue",
    "create_http_record_sender",
    "group_records_by_type",
]
