"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/10
@description: transport 层上传结果类型定义
"""

from dataclasses import dataclass, field
from typing import List, Optional, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


@dataclass
class UploadFailure:
    """上传失败的诊断信息。"""

    message: str
    records: List[Record] = field(default_factory=list)
    status: Optional[int] = None
    code: Optional[Union[int, str]] = None
    exception: Optional[BaseException] = None


@dataclass
class UploadResult:
    """单次上传或分发汇总结果。"""

    # confirmed_records 表示本轮已终结、不再重试的记录，包括成功记录和已报告的确定性失败记录。
    confirmed_records: List[Record] = field(default_factory=list)
    remaining_records: List[Record] = field(default_factory=list)
    failures: List[UploadFailure] = field(default_factory=list)


__all__ = [
    "UploadFailure",
    "UploadResult",
]
