"""
@author: caddiesnew
@file: buffer.py
@time: 2026/4/17
@description: 基于 list + set 去重索引的 RecordBuffer
"""

from typing import List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


class RecordBuffer:
    """
    Record 缓冲区，内部用 list 存储 Record、set[int] 做 num 去重索引。

    调用方需要在外部加锁（如 threading.Condition）。
    """

    __slots__ = ("_records", "_record_num_index")

    def __init__(self) -> None:
        self._records: List[Record] = []
        self._record_num_index: set[int] = set()

    def __len__(self) -> int:
        return len(self._records)

    def __bool__(self) -> bool:
        return bool(self._records)

    # ── 写入 ──

    def extend(self, records: List[Record]) -> int:
        """追加 records，自动按 num 去重。返回实际入队的数量。"""
        accepted_records = [record for record in records if self._try_enqueue(record)]
        self._records.extend(accepted_records)
        return len(accepted_records)

    def prepend(self, records: List[Record]) -> int:
        """回滚到头部，自动按 num 去重。返回实际入队的数量。"""
        accepted_records = [record for record in records if self._try_enqueue(record)]
        if accepted_records:
            self._records[:0] = accepted_records
        return len(accepted_records)

    def _try_enqueue(self, record: Record) -> bool:
        """根据 record num 去重，未存在则注册并返回 True。"""
        if record.num in self._record_num_index:
            return False
        self._record_num_index.add(record.num)
        return True

    # ── 读取 ──

    def drain(self) -> List[Record]:
        """取出全部 records 并清空缓冲区（含索引）。"""
        pending_records = self._records[:]
        self._records.clear()
        self._record_num_index.clear()
        return pending_records


__all__ = [
    "RecordBuffer",
]
