"""
@author: nexisato
@file: __init__.py
@time: 2026/10/23
@description: sync 范围过滤器，按 step 与时间裁剪本地回放记录
"""

from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import SyncRange

__all__ = ["SyncRangeFilter"]


class SyncRangeFilter:
    """按 step 范围与时间范围过滤回放记录，任一维度未设置表示该维度不限制。

    step 范围仅对 scalar/media 记录生效，时间范围对 scalar/media/log 记录生效，
    两个维度同时设置时取交集，上下界均为闭区间。
    """

    __slots__ = ("_start_step", "_end_step", "_start_ns", "_end_ns")

    def __init__(
        self,
        *,
        start_step: Optional[int] = None,
        end_step: Optional[int] = None,
        start_ns: Optional[int] = None,
        end_ns: Optional[int] = None,
    ):
        self._start_step = start_step
        self._end_step = end_step
        self._start_ns = start_ns
        self._end_ns = end_ns

    @classmethod
    def from_proto(cls, sync_range: Optional[SyncRange]) -> "SyncRangeFilter":
        """
        从协议对象构建过滤器，None 表示全量同步。

        时间戳为 0（未设置或 epoch）时视为不限制，避免把未赋值字段当作 1970 年过滤。
        """
        if sync_range is None:
            return cls()
        start_step = sync_range.start_step if sync_range.HasField("start_step") else None
        end_step = sync_range.end_step if sync_range.HasField("end_step") else None
        start_ns = sync_range.start_time.ToNanoseconds() if sync_range.HasField("start_time") else None
        end_ns = sync_range.end_time.ToNanoseconds() if sync_range.HasField("end_time") else None
        return cls(
            start_step=start_step,
            end_step=end_step,
            start_ns=start_ns if start_ns else None,
            end_ns=end_ns if end_ns else None,
        )

    @property
    def enabled(self) -> bool:
        """是否启用了任意维度的范围过滤。"""
        return any(value is not None for value in (self._start_step, self._end_step, self._start_ns, self._end_ns))

    @property
    def has_time(self) -> bool:
        """是否设置了时间范围。终端日志只受时间范围约束，未设置时不参与范围回放。"""
        return self._start_ns is not None or self._end_ns is not None

    def match(self, *, step: Optional[int] = None, timestamp: Optional[Timestamp] = None) -> bool:
        """
        判断记录是否落在范围内。

        :param step: 记录 step，None 表示该记录没有 step（如终端日志），不参与 step 过滤
        :param timestamp: 记录时间，None 表示不参与时间过滤
        :return: True 表示保留，False 表示跳过
        """
        if step is not None:
            if self._start_step is not None and step < self._start_step:
                return False
            if self._end_step is not None and step > self._end_step:
                return False
        if timestamp is not None:
            ns = timestamp.ToNanoseconds()
            if self._start_ns is not None and ns < self._start_ns:
                return False
            if self._end_ns is not None and ns > self._end_ns:
                return False
        return True
