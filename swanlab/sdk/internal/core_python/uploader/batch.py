"""
@author: caddiesnew
@file: batch.py
@time: 2026/3/19
@description: Record 批处理辅助函数
"""

from typing import Iterator, Sequence, Tuple, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

RecordLike = Union[Record, bytes]


def load_record(record_like: RecordLike) -> Record:
    """统一接收 protobuf Record 或其序列化 bytes。"""
    if isinstance(record_like, Record):
        return record_like
    if isinstance(record_like, bytes):
        record = Record()
        record.ParseFromString(record_like)
        return record
    raise TypeError(f"Unsupported record type: {type(record_like).__name__}")


def _generate_chunks(records: Sequence[RecordLike], per_request_len: int) -> Iterator[Tuple[Sequence[RecordLike], int]]:
    """
    按固定大小对 Record 序列做切片。
    yield: (分片数据, 分片长度)
    """
    if per_request_len == -1:
        yield records, len(records)
        return
    if per_request_len <= 0:
        raise ValueError(f"per_request_len must be -1 or > 0, got {per_request_len}")

    total = len(records)
    if total <= per_request_len:
        yield records, total
        return

    for index in range(0, total, per_request_len):
        chunk = records[index : index + per_request_len]
        yield chunk, len(chunk)


__all__ = [
    "RecordLike",
    "load_record",
]
