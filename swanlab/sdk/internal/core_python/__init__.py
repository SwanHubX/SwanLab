"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: SwanLab Core 的 Python sidecar 适配层。

实现 CoreProtocol，当前为纯 Python 实现：
1. 本地落盘仍在 Python 内完成；
2. 上传线程负责缓冲、聚合和重试；
3. 与后端的交互未来统一下沉到 swanlab-core（Go sidecar）。

在 swanlab-core 尚未接入前，上传 transport 使用占位实现，
这样 BackgroundConsumer 等调用方可以先稳定对接 CoreProtocol，
后续替换为 Go Core 时无需修改业务入口。
"""

from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core import CoreProtocol
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.uploader import Uploader
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.utils.helper.env import DEBUG

__all__ = ["CorePython"]


class CorePython(CoreProtocol):
    """
    CoreProtocol 的 Python 实现。
    由 SwanLabRun 在初始化时构造并注入给 BackgroundConsumer。
    """

    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        self._store: Optional[DataStoreWriter] = None
        self._uploader: Optional[Uploader] = None

    def startup(self, cloud: bool, persistence: bool) -> None:
        if self._store is not None or self._uploader is not None:
            raise RuntimeError("CorePython has already been started.")
        if persistence:
            self._store = DataStoreWriter()
            self._store.open(str(self._ctx.run_file))
        if cloud:
            self._uploader = Uploader()

    def handle_records(self, records: List[Record]) -> None:
        if self._store is None and self._uploader is None:
            console.warning("CorePython is not started, skipping record handling.")
            return
        if self._store is not None:
            for record in records:
                self._store.write(record.SerializeToString())
                if DEBUG:
                    console.debug("Write record:", record.WhichOneof("record_type"))
        if self._uploader is not None:
            self._uploader.put(records)

    def shutdown(self) -> None:
        if self._uploader is not None:
            self._uploader.finish()
            self._uploader = None
        if self._store is not None:
            self._store.close()
            self._store = None
