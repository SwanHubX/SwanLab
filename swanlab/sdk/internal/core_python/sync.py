"""
@author: cunyue
@file: sync.py
@time: 2026/5/14 14:21
@description: sync 协议实现
"""

from typing import Optional

from swanlab.exceptions import DataStoreError
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest, DeliverSyncStartResponse
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg import executor
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.core_python.transport import Transport
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.protocol.core import CoreSyncProtocol

__all__ = ["CoreSyncPython"]


class CoreSyncPython(CoreSyncProtocol):
    def __init__(self):
        super().__init__()
        self._run_ctx: Optional[CoreContext] = None
        self._reader = DataStoreReader()
        self._read_executor = executor.EventLoopExecutor("SwanLab·Sync·Reader")
        self._transport: Optional[Transport] = None

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        self._ctx = CoreContext.from_proto(start_request.core_settings)
        # 安全地开启sync，DataStoreError 通常是run头文件损坏，属于已知错误
        with safe.block(message="Failed to open sync store with unexpected error", write_to_tty=False):
            try:
                self._reader.open(self._ctx.run_file)
                self._read_executor.start(self.read())
                return DeliverSyncStartResponse(success=True, message="success")
            except DataStoreError as e:
                return DeliverSyncStartResponse(success=False, message=str(e))
        return DeliverSyncStartResponse(success=False, message="unknown error")

    async def read(self):
        """
        异步读取本地文件，存储在内存中
        """
        self._transport = Transport(self._ctx)
        for record_bytes in self._reader:
            record = Record()
            record.ParseFromString(record_bytes)
            self._transport.put(record)
            # TODO 处理traker记录

    def deliver_sync_finish(self):
        # 1. 等待读取线程结束
        self._read_executor.wait()
        self._reader.close()
        # 2. 云端重新开启

    async def upload(self):
        """
        异步上传处理后的记录
        """
        ...
