"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传抽象（HTTP 传输层）
"""

from typing import Optional, Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console


class HttpRecordSender:
    """PLACEHOLDER: HTTP 传输层上传 sender，请求映射待实现。"""

    def upload(self, record_type: str, records: Sequence[Record]) -> None:
        """通用上传入口，按 record_type 路由到对应 upload_{kind} 方法。"""
        if len(records) == 0:
            return
        fn = getattr(self, f"upload_{record_type}", None)
        if fn is not None:
            fn(records)

    def upload_run(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_run (request mapping pending).")

    def upload_finish(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_finish (request mapping pending).")

    def upload_column(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_column (request mapping pending).")

    def upload_metric(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_metric (request mapping pending).")

    def upload_config(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_config (request mapping pending).")

    def upload_console(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_console (request mapping pending).")

    def upload_metadata(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_metadata (request mapping pending).")

    def upload_requirements(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_requirements (request mapping pending).")

    def upload_conda(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_conda (request mapping pending).")

    def close(self) -> None:
        pass


# ── 模块级单例 ──

_default_sender: Optional[HttpRecordSender] = None


def create_record_sender() -> HttpRecordSender:
    """获取或创建全局 Record sender 单例。"""
    global _default_sender
    if _default_sender is None:
        _default_sender = HttpRecordSender()
    return _default_sender


def reset_record_sender() -> None:
    """重置全局 Record sender 单例。"""
    global _default_sender
    if _default_sender is not None:
        _default_sender.close()
        _default_sender = None


__all__ = [
    "HttpRecordSender",
    "create_record_sender",
    "reset_record_sender",
]
