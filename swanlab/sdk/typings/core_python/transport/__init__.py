"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/10
@description: transport 层上传结果类型定义
"""

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

# 失败来源：说明"哪一层失败了"，用于诊断，不参与重试决策
FailureSource = Literal["swanlab_api", "object_storage", "local_data", "worker_internal"]
# - swanlab_api：SwanLab 后端接口，如 house/metrics、columns、profile、prepare / complete
# - object_storage：对象存储字节传输，如媒体 PUT、save 文件与分片 PUT
# - local_data：本地数据不可用，如文件缺失/不可读/超限/解析失败
# - worker_internal：dispatch 未预期异常，无法确认请求是否完成

# 失败阶段：说明"流程走到哪一步失败了"
FailureStage = Literal["request_address", "transfer", "submit_metrics", "complete"]
# - request_address：申请上传地址（save prepare、预签名数量校验）
# - transfer：对象存储字节传输（媒体与 save 的 PUT，含分片；媒体预签名请求失败也归入此阶段）
# - submit_metrics：提交指标与配置（house/metrics、columns、profile）
# - complete：确认上传完成（complete / complete-multipart）


@dataclass
class UploadFailure:
    """上传失败的诊断信息。"""

    source: FailureSource
    message: str
    records: List[Record] = field(default_factory=list)
    stage: Optional[FailureStage] = None
    status: Optional[int] = None
    code: Optional[Union[int, str]] = None
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

    # confirmed_records 表示本轮已终结、不再重试的记录，包括成功记录和已报告的确定性失败记录。
    confirmed_records: List[Record] = field(default_factory=list)
    remaining_records: List[Record] = field(default_factory=list)
    failures: List[UploadFailure] = field(default_factory=list)


__all__ = [
    "FailureSource",
    "FailureStage",
    "UploadFailure",
    "UploadResult",
]
