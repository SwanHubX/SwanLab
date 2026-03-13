"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: SwanLab Core Python 版本，封装SwanLab云端版核心业务，包括：
1. 提供http客户端，用于与SwanLab云端API进行交互。
2. 提供rpc封装函数，以rpc方式调用SwanLab云端API。
3. 提供上传线程，在另一个线程执行上传任务。
...

实现 CoreProtocol，当前为纯 Python 实现。
未来由 swanlab-core（Go 二进制）替代时，此模块整体被替换，
BackgroundConsumer 等调用方无需修改。
"""

from typing import List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core import CoreProtocol
from swanlab.sdk.internal.core_python.store import DataStoreWriter
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
        self._store: DataStoreWriter | None = None

    def startup(self, cloud: bool, persistence: bool) -> None:
        if self._store is not None:
            raise RuntimeError("CorePython has already been started.")
        if persistence:
            self._store = DataStoreWriter()
            self._store.open(str(self._ctx.run_file))

    def handle_records(self, records: List[Record]) -> None:
        if self._store is None:
            console.warning("CorePython is not started, skipping record handling.")
            return
        for record in records:
            self._store.write(record.SerializeToString())
            if DEBUG:
                console.debug("[CORE-PY] Write record:", record.WhichOneof("record_type"))

    def shutdown(self) -> None:
        if self._store is None:
            return
        self._store.close()
        self._store = None
