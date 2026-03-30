#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@author: caddiesnew
@file: helper.py
@time: 2026/3/30 15:33
@description: uploader 辅助函数
"""

from collections import OrderedDict
from typing import Dict, Iterator, List, Sequence, Tuple

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


def generate_chunks(records: Sequence[Record], per_request_len: int) -> Iterator[Tuple[Sequence[Record], int]]:
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


def group_records_by_type(records: Sequence[Record]) -> Dict[str, List[Record]]:
    grouped: Dict[str, List[Record]] = OrderedDict()
    for record in records:
        if not isinstance(record, Record):
            raise TypeError(f"group_records_by_type only accepts Record instances, got {type(record).__name__}")
        record_type = record.WhichOneof("record_type")
        if record_type is None:
            raise ValueError("Record has no active record_type and cannot be uploaded.")
        grouped.setdefault(record_type, []).append(record)
    return grouped
