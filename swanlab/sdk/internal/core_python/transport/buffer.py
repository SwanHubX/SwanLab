"""
@author: caddiesnew
@file: buffer.py
@time: 2026/4/17
@description: 按编号去重的 RecordBuffer
"""

from typing import Iterable, List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


class RecordBuffer:
    """
    Record 缓冲区，按插入顺序存储记录，同编号保留最后写入的内容。

    调用方需要在外部加锁（如 threading.Condition）。
    """

    __slots__ = ("_records",)

    def __init__(self) -> None:
        self._records: dict[int, Record] = {}

    def __len__(self) -> int:
        return len(self._records)

    def __bool__(self) -> bool:
        return bool(self._records)

    # ── 写入 ──

    def extend(self, records: Iterable[Record]) -> int:
        """覆盖同编号记录，保留原位置。返回新增记录数。"""
        previous_size = len(self._records)
        for record in records:
            self._records[record.num] = record
        return len(self._records) - previous_size

    def prepend(self, records: List[Record]) -> int:
        """回滚到头部，保留缓冲中的同编号记录。返回新增记录数。"""
        accepted_records: dict[int, Record] = {}
        for record in records:
            if record.num not in self._records:
                accepted_records[record.num] = record
        accepted = len(accepted_records)
        if accepted:
            accepted_records.update(self._records)
            self._records = accepted_records
        return accepted

    # ── 读取 ──

    def drain(self) -> List[Record]:
        """取出全部 records 并清空缓冲区（含索引）。"""
        pending_records = list(self._records.values())
        self._records.clear()
        return pending_records


__all__ = [
    "RecordBuffer",
]
