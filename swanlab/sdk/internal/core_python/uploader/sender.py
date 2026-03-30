"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/29
@description: HTTP uploader grouping and sender abstractions
"""

from collections import defaultdict
from typing import DefaultDict, Protocol

from swanlab.sdk.internal.bus.events import FlushPayload, MetricsUploadEvent, MetricsUploadPayload
from swanlab.sdk.internal.pkg import console


def group_records_by_type(records: FlushPayload) -> MetricsUploadPayload:
    """Group protobuf records by their oneof record_type."""
    grouped: DefaultDict[str, list] = defaultdict(list)
    for record in records:
        record_type = record.WhichOneof("record_type")
        if record_type is None:
            raise ValueError("Record must contain a record_type before upload.")
        grouped[record_type].append(record)
    return dict(grouped)


class HttpRecordSender(Protocol):
    """HTTP sender used by the batch uploader."""

    def send(self, event: MetricsUploadEvent) -> None: ...

    def close(self) -> None: ...


class NoopHttpRecordSender:
    """
    Placeholder HTTP sender.

    The SDK already persists protobuf records locally before enqueueing them here,
    so a temporary no-op sender will not affect data durability.
    """

    def __init__(self) -> None:
        self._announced = False

    def send(self, event: MetricsUploadEvent) -> None:
        if self._announced:
            return
        console.debug(
            "CorePython HTTP uploader is buffering protobuf records, but the remote sender is not implemented yet."
        )
        self._announced = True

    def close(self) -> None:
        pass


def create_http_record_sender() -> HttpRecordSender:
    """Create the HTTP sender used by CorePython cloud uploads."""
    return NoopHttpRecordSender()


__all__ = [
    "HttpRecordSender",
    "NoopHttpRecordSender",
    "create_http_record_sender",
    "group_records_by_type",
]
