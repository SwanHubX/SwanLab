"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/10
@description: transport 层共享结果类型：一次上传尝试（sender）或一次汇总（dispatch）的结构化结果

UploadResult / UploadFailure 是 worker、dispatch、sender 之间的确认契约：
- confirmed / remaining 以 Record 对象身份区分，同一记录不得同时出现在两侧
- 失败与未确认记录关联；并发文件上传允许产生多个失败
- exception 仅供内部诊断（日志文件/调试），不得整体打印到终端（可能包含 URL 等信息）
"""

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

# 失败来源：SwanLab API、对象存储、本地数据处理、worker 内部错误
FailureSource = Literal["swanlab_api", "object_storage", "local_data", "worker_internal"]

# 操作阶段：请求上传地址、传输内容、提交指标、complete 确认
FailureStage = Literal["request_address", "transfer", "submit_metrics", "complete"]


@dataclass
class UploadFailure:
    """一次上传失败的结构化描述，字段按实际来源填写。

    message 与 summary 只承载安全诊断信息（来源、阶段、HTTP 状态、业务错误码），
    不包含预签名 URL、Authorization、完整响应或原始异常文本。
    """

    source: FailureSource
    message: str
    # 关联的未确认记录；批量级失败关联整批，文件级失败关联单条
    records: List[Record] = field(default_factory=list)
    stage: Optional[FailureStage] = None
    status: Optional[int] = None
    code: Optional[Union[int, str]] = None
    # 可解析的 Retry-After（秒），由重试调度填充
    retry_after: Optional[float] = None
    # 原始异常，仅供内部诊断，不写入 proto、不直接整体打印
    exception: Optional[BaseException] = None

    @property
    def summary(self) -> str:
        """安全的一行失败摘要，用于节流告警与停滞诊断。"""
        parts = [self.source, self.stage or "-"]
        if self.status is not None:
            parts.append(f"http {self.status}")
        if self.code is not None:
            parts.append(f"code={self.code}")
        return " ".join(parts)


@dataclass
class UploadResult:
    """一次 sender 尝试或一次 dispatch 汇总的结果。

    契约（plan §4.1）：
    1. 输入中的有效记录均可归入 confirmed 或 remaining，不凭空创建、不因失败消失
    2. 同一记录不得同时处于 confirmed 与 remaining
    3. 无失败返回不代表成功，是否消费由 confirmed 明确决定
    """

    confirmed_records: List[Record] = field(default_factory=list)
    remaining_records: List[Record] = field(default_factory=list)
    failures: List[UploadFailure] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """本轮任务是否全部确认。"""
        return not self.remaining_records and not self.failures


__all__ = [
    "FailureSource",
    "FailureStage",
    "UploadFailure",
    "UploadResult",
]
