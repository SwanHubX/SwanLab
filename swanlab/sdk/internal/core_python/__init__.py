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
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRequest, FinishResponse, StartRequest, StartResponse
from swanlab.sdk.internal.context import CallbackManager, RunContext
from swanlab.sdk.protocol import CoreProtocol

__all__ = ["CorePython"]


class CorePython(CoreProtocol):
    """
    CoreProtocol 的 Python 实现。
    由 Run 在初始化时构造并注入给 BackgroundConsumer。
    """

    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        self._callbacker: CallbackManager = ctx.callbacker

    def start(self, start_request: StartRequest) -> StartResponse: ...

    def publish(self, records: List[Record]) -> None:
        pass

    def fork(self) -> "CorePython":
        raise NotImplementedError(
            "CorePython.fork() is not implemented. Please waiting for go version, while you should not reach here?"
        )

    def finish(self, finish_request: FinishRequest) -> FinishResponse: ...
