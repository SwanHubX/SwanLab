"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 上传结果汇总
"""

from collections import OrderedDict
from typing import TYPE_CHECKING, Dict, Iterator, List, Optional, Sequence, Tuple

from swanlab.exceptions import NonRetryableUploadError
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.core_python.transport import FailureSource, UploadFailure, UploadResult

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


def generate_chunks(records: Sequence[Record], batch_size: int) -> Iterator[Tuple[Sequence[Record], int]]:
    """
    按固定大小对 Record 序列做切片。
    yield: (分片数据, 分片长度)
    """
    if not records:
        return
    if batch_size == -1:
        yield records, len(records)
        return
    if batch_size <= 0:
        raise ValueError(f"max_records_per_request must be -1 or > 0, got {batch_size}")

    total = len(records)
    for index in range(0, total, batch_size):
        chunk = records[index : index + batch_size]
        yield chunk, len(chunk)


def group_records_by_type(records: Sequence[Record]) -> Dict[str, List[Record]]:
    """按 record_type 对 Record 序列分组，保持插入顺序。"""
    records_by_type: Dict[str, List[Record]] = OrderedDict()
    for record in records:
        if not isinstance(record, Record):
            raise TypeError(f"Expected Record instance, got {type(record).__name__}")
        record_type = record.WhichOneof("record_type")
        if record_type is None:
            raise ValueError("Failed to determine record type, maybe you are using an incompatible proto version?")
        records_by_type.setdefault(record_type, []).append(record)
    return records_by_type


class Dispatch:
    """按 record_type 聚合记录并分批上传，汇总上传结果与记录数统计。"""

    _RECORD_TYPES = frozenset(f.name for f in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)

    def __init__(
        self,
        batch_size: int,
        sender: Optional[HttpRecordSender] = None,
        tracker: Optional["UploadTracker"] = None,
    ):
        self._batch_size = batch_size
        self._sender = sender
        self._tracker = tracker

    def __call__(self, records: List[Record]) -> UploadResult:
        """按 record_type 聚合并上传，返回汇总结果。"""
        try:
            records_by_type = group_records_by_type(records)
        except (ValueError, TypeError) as e:
            # 本地 record 格式异常丢弃，避免被 worker 重试
            return self._settle(records, str(e))
        confirmed: List[Record] = []
        failures: List[UploadFailure] = []
        record_types = list(records_by_type.keys())
        for index, record_type in enumerate(record_types):
            result = self._upload_record_type(record_type, records_by_type[record_type])
            confirmed.extend(result.confirmed_records)
            failures.extend(result.failures)
            if result.remaining_records:
                unprocessed_records = [
                    record
                    for next_record_type in record_types[index + 1 :]
                    for record in records_by_type[next_record_type]
                ]
                return UploadResult(
                    confirmed_records=confirmed,
                    remaining_records=result.remaining_records + unprocessed_records,
                    failures=failures,
                )
        return UploadResult(confirmed_records=confirmed, failures=failures)

    def _upload_record_type(self, record_type: str, records: List[Record]) -> UploadResult:
        """上传单个 record_type 分组。"""
        if record_type not in self._RECORD_TYPES:
            console.warning(f"No handler for record kind={record_type!r}, skipping {len(records)} record(s).")
            return self._settle(records, f"unknown record type {record_type!r}")
        confirmed: List[Record] = []
        failures: List[UploadFailure] = []
        offset = 0
        for chunk, chunk_size in generate_chunks(records, self._batch_size):
            result = self._upload_chunk(record_type, chunk)
            if len(result.confirmed_records) + len(result.remaining_records) != chunk_size:
                console.error(f"Upload result did not account for all {chunk_size} record(s) of {record_type!r}.")
                failures.append(
                    UploadFailure(
                        source="worker_internal",
                        message=f"upload result did not account for all {chunk_size} record(s)",
                        records=list(chunk),
                    )
                )
                return UploadResult(
                    confirmed_records=confirmed, remaining_records=list(records[offset:]), failures=failures
                )
            self._advance_records(len(result.confirmed_records))
            confirmed.extend(result.confirmed_records)
            failures.extend(result.failures)
            if result.remaining_records:
                return UploadResult(
                    confirmed_records=confirmed,
                    remaining_records=result.remaining_records + records[offset + chunk_size :],
                    failures=failures,
                )
            offset += chunk_size
        return UploadResult(confirmed_records=confirmed, failures=failures)

    def _upload_chunk(self, record_type: str, chunk: Sequence[Record]) -> UploadResult:
        """上传单个分片，未预期异常转为当前分片失败结果。"""
        if self._sender is None:
            raise RuntimeError("sender not set")
        if record_type == "log":
            console.debug(
                f"Dispatching log records: count={len(chunk)}, nums={chunk[0].num}..{chunk[-1].num}",
                write_to_tty=False,
            )
        try:
            result = self._sender.upload(record_type, chunk)
        except NonRetryableUploadError as e:
            # 确定性失败：终结当前分片，避免阻塞同批后续 record 的上传
            return UploadResult(
                confirmed_records=list(chunk),
                failures=[UploadFailure(source="swanlab_api", message=str(e), records=list(chunk), exception=e)],
            )
        except Exception as e:
            source: FailureSource = "local_data" if isinstance(e, OSError) else "worker_internal"
            return UploadResult(
                confirmed_records=list(chunk) if source == "local_data" else [],
                remaining_records=[] if source == "local_data" else list(chunk),
                failures=[
                    UploadFailure(
                        source=source,
                        message=f"unexpected {type(e).__name__} while uploading {record_type!r} chunk",
                        records=list(chunk),
                        exception=e,
                    )
                ],
            )
        if record_type == "log":
            console.debug(
                f"Log sender completed: count={len(chunk)}, nums={chunk[0].num}..{chunk[-1].num}",
                write_to_tty=False,
            )
        return result

    def _settle(self, records: Sequence[Record], message: str) -> UploadResult:
        """本地数据无法处理：终结整批记录并推进统计。"""
        self._advance_records(len(records))
        return UploadResult(
            confirmed_records=list(records),
            failures=[UploadFailure(source="local_data", message=message, records=list(records))],
        )

    def _advance_records(self, count: int) -> None:
        """根据确认记录数更新 tracker 统计。"""
        if self._tracker is not None and count > 0:
            self._tracker.advance_records(count)


__all__ = [
    "Dispatch",
    "generate_chunks",
    "group_records_by_type",
]
