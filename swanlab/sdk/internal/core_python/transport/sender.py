"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传抽象（HTTP 传输层）
"""

from typing import Sequence

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

    def upload_start(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_start (request mapping pending).")

    def upload_run(self, records: Sequence[Record]) -> None:
        self.upload_start(records)

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


__all__ = [
    "HttpRecordSender",
]
