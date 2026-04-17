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

from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRequest, FinishResponse, StartRequest, StartResponse
from swanlab.sdk.internal.context import CallbackManager, RunContext
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.transport import Transport, reset_record_sender
from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.protocol import CoreProtocol

__all__ = ["CorePython"]


class CorePython(CoreProtocol):
    """
    CoreProtocol 的 Python 实现。
    由 Run 在初始化时构造并注入给 BackgroundConsumer。
    """

    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        self._store: Optional[DataStoreWriter] = None
        self._transport: Optional[Transport] = None
        self._callbacker: CallbackManager = ctx.callbacker
        self._mode = ctx.config.settings.mode

    def start(self, start_request: StartRequest) -> StartResponse:
        if self._store is not None or self._transport is not None:
            raise RuntimeError("CorePython has already been started.")
        if self._mode != "disabled":
            self._store = DataStoreWriter()
            self._store.open(str(self._ctx.run_file))
        if self._mode == "cloud":
            self._transport = Transport()
        return StartResponse(success=True, color="#ffffff")

    def publish(self, records: List[Record]) -> None:
        if self._store is None and self._transport is None:
            console.warning("CorePython is not started, skipping record handling.")
            return
        if self._store is not None:
            for record in records:
                self._store.write(record.SerializeToString())
                if helper.DEBUG:
                    console.debug("Write record:", record.WhichOneof("record_type"))
        if self._transport is not None:
            self._transport.put(records)

    def fork(self) -> "CorePython":
        raise NotImplementedError(
            "CorePython.fork() is not implemented. Please waiting for go version, while you should not reach here?"
        )

    def finish(self, finish_request: FinishRequest) -> FinishResponse:
        if self._transport is not None:
            self._transport.finish()
            self._transport = None
        reset_record_sender()
        if self._store is not None:
            self._store.close()
            self._store = None

        return FinishResponse(success=True, message="I'm not ready.")
