"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:30
@description: SwanLab 核心模块
"""

from typing import TYPE_CHECKING, List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, FinishResponse, StartRecord, StartResponse
from swanlab.sdk.protocol import CoreEnum, CoreProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


class NullCore(CoreProtocol):
    """空 Core，所有方法为 no-op"""

    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        return StartResponse(success=True, message="I'm a teapot.", run=start_record)

    def publish(self, records: List[Record]) -> None: ...

    def fork(self) -> "NullCore":
        raise NotImplementedError("NullCore does not support fork, you should not reach here?")

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        return FinishResponse(success=True, message="I'm a teapot.")


# TODO: 未来实现core以后，python版本依旧会有一段时间的同时存在时间。后续实现一种机制，选择不同的core实现
core_enum: CoreEnum = CoreEnum.CORE_PYTHON


def create_core(ctx: "RunContext") -> CoreProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if ctx.config.settings.mode == "disabled":
        return NullCore(ctx)

    if core_enum == CoreEnum.CORE_PYTHON:
        from swanlab.sdk.internal.core_python import CorePython

        return CorePython(ctx)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError(f"CoreEnum {core_enum} is not supported yet.")
