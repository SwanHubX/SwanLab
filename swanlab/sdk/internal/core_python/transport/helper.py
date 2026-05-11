#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@author: caddiesnew
@file: helper.py
@time: 2026/3/30 15:33
@description: transport 辅助函数
"""

import hashlib
from collections import OrderedDict
from typing import Dict, Iterator, List, Sequence, Tuple, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


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
    if total <= batch_size:
        yield records, total
        return

    for index in range(0, total, batch_size):
        chunk = records[index : index + batch_size]
        yield chunk, len(chunk)


def group_records_by_type(records: Sequence[Record]) -> Dict[str, List[Record]]:
    """
    使用 WhichOneof("record_type") 按 record 类型分组。

    record_type = record.WhichOneof("record_type") 返回 oneof 中实际设置的字段名，
    如 "start"、"metric" 等，与 proto 定义一一对应。

    返回 OrderedDict 保持插入顺序。
    """
    records_by_type: Dict[str, List[Record]] = OrderedDict()
    for record in records:
        if not isinstance(record, Record):
            raise TypeError(f"Expected Record instance, got {type(record).__name__}")
        record_type = record.WhichOneof("record_type")
        if record_type is None:
            raise ValueError("Record.record_type is not set")
        records_by_type.setdefault(record_type, []).append(record)
    return records_by_type


def compute_md5(source: Union[str, bytes], chunk_size: int = 8 * 1024 * 1024) -> str:
    """计算 MD5。接受文件路径（str）或已读取的字节数据（bytes）。"""
    md5 = hashlib.md5()
    if isinstance(source, bytes):
        md5.update(source)
    else:
        with open(source, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5.update(chunk)
    return md5.hexdigest()
