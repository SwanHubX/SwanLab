"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/10
@description: transport 层上传结果类型定义
"""

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

FailureSource = Literal["swanlab_api", "object_storage", "local_data", "worker_internal"]

FailureStage = Literal["request_address", "transfer", "submit_metrics", "complete"]


@dataclass
class UploadFailure:
    """上传失败的诊断信息。"""

    source: FailureSource
    message: str
    records: List[Record] = field(default_factory=list)
    stage: Optional[FailureStage] = None
    status: Optional[int] = None
    code: Optional[Union[int, str]] = None
    retry_after: Optional[float] = None
    exception: Optional[BaseException] = None

    @property
    def summary(self) -> str:
        """单行失败摘要。"""
        parts = [self.source, self.stage or "-"]
        if self.status is not None:
            parts.append(f"http {self.status}")
        if self.code is not None:
            parts.append(f"code={self.code}")
        return " ".join(parts)


@dataclass
class UploadResult:
    """单次上传或分发汇总结果。"""

    confirmed_records: List[Record] = field(default_factory=list)
    remaining_records: List[Record] = field(default_factory=list)
    failures: List[UploadFailure] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """本轮任务是否完成。"""
        return not self.remaining_records and not self.failures


__all__ = [
    "FailureSource",
    "FailureStage",
    "UploadFailure",
    "UploadResult",
]
