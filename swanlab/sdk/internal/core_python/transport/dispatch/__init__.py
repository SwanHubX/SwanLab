"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 上传结果汇总
"""

from collections import OrderedDict
from typing import TYPE_CHECKING, Dict, Iterator, List, Optional, Sequence, Tuple

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

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


def group_records_by_type(records: Sequence[Record]) -> Tuple[Dict[str, List[Record]], List[Record]]:
    """按 record_type 对 Record 序列分组，保持插入顺序。

    返回 (分组结果, 无法识别类型的记录)。后者已终止处理，由调用方计入 confirmed，
    不影响同批其他记录的上传。
    """
    records_by_type: Dict[str, List[Record]] = OrderedDict()
    invalid_records: List[Record] = []
    for record in records:
        if not isinstance(record, Record):
            raise TypeError(f"Expected Record instance, got {type(record).__name__}")
        record_type = record.WhichOneof("record_type")
        if record_type is None:
            invalid_records.append(record)
            continue
        records_by_type.setdefault(record_type, []).append(record)
    return records_by_type, invalid_records


class Dispatch:
    """按 record_type 聚合记录并分批上传，汇总上传结果与记录数统计。"""

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
        records_by_type, invalid_records = group_records_by_type(records)
        confirmed: List[Record] = list(invalid_records)
        failures: List[UploadFailure] = []
        if invalid_records:
            # 无法识别类型的记录已终止处理：跳过上传但计入 confirmed
            console.warning(f"Failed to determine record type, skipping {len(invalid_records)} record(s).")
            self._advance_records(len(invalid_records))
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
        confirmed: List[Record] = []
        failures: List[UploadFailure] = []
        offset = 0
        for chunk, chunk_size in generate_chunks(records, self._batch_size):
            result = self._upload_chunk(record_type, chunk)
            if len(result.confirmed_records) + len(result.remaining_records) != chunk_size:
                console.error(f"Upload result did not account for all {chunk_size} record(s) of {record_type!r}.")
                failures.append(
                    UploadFailure(
                        message=f"upload result did not account for all {chunk_size} record(s)",
                        records=list(chunk),
                    )
                )
                return UploadResult(
                    confirmed_records=confirmed, remaining_records=list(records[offset:]), failures=failures
                )
            # uploaded_records 表示“已终止处理 Record 数”，confirmed 均计入
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
        """上传单个分片，兜底 sender 未处理的未预期异常。

        正常的业务/网络错误由 sender 转为 UploadResult，不会进入这里的 except。
        兜底策略：
        - 不可重试的 ApiError：确定性业务拒绝，终结当前分片；
        - 非 ApiError 的 OSError：本地数据问题，重试不会恢复，终结；
        - 其它未知异常：fail-open 保留重试，避免整片记录被误丢。
        """
        if self._sender is None:
            raise RuntimeError("sender not set")
        if record_type == "log":
            console.debug(
                f"Dispatching log records: count={len(chunk)}, nums={chunk[0].num}..{chunk[-1].num}",
                write_to_tty=False,
            )
        try:
            result = self._sender.upload(record_type, chunk)
        except Exception as e:
            # sender.upload 会捕获所有的 requests.RequestException，本本部分仅做兜底处理
            # upload 的捕获行为在此处通过注释做隐式约定：如果没有捕获 RequestException ，本部分会直接丢掉这部分网络错误，重试行为实效
            # 1. 后端明确拒绝（通常是 4xx）：终结当前分片，避免阻塞同批后续 record 的上传
            if isinstance(e, ApiError) and not e.retryable:
                return UploadResult(
                    confirmed_records=list(chunk),
                    failures=[UploadFailure(message=str(e), records=list(chunk), exception=e)],
                )
            # 2. 根据异常类型决定终结还是重试
            local_failure = isinstance(e, OSError) and not isinstance(e, ApiError)
            if local_failure:
                # 2.1. 非 HTTP 的 OSError 视为本地数据问题：重试不会恢复，终结当前分片
                # 注意：requests.RequestException 也继承 OSError；当前 sender.upload 已把它转为 UploadResult，因此这里只应收到非 HTTP 的本地 OSError。
                confirmed_records = list(chunk)
                remaining_records = []
            else:
                # 2.2. 其它未知异常：fail-open 保留重试，避免整片记录被误丢，代价是确定性的未知异常可能长时间重试。
                confirmed_records = []
                remaining_records = list(chunk)

            failure = UploadFailure(
                message=f"unexpected {type(e).__name__} while uploading {record_type!r} chunk",
                records=list(chunk),
                exception=e,
            )
            return UploadResult(
                confirmed_records=confirmed_records, remaining_records=remaining_records, failures=[failure]
            )
        if record_type == "log":
            console.debug(
                f"Log sender completed: count={len(chunk)}, nums={chunk[0].num}..{chunk[-1].num}",
                write_to_tty=False,
            )
        return result

    def _advance_records(self, count: int) -> None:
        """根据确认记录数更新 tracker 统计。"""
        if self._tracker is not None and count > 0:
            self._tracker.advance_records(count)


__all__ = [
    "Dispatch",
    "generate_chunks",
    "group_records_by_type",
]
