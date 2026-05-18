"""
@author: cunyue
@file: sync.py
@time: 2026/5/14 14:21
@description: sync 协议实现
"""

from typing import Callable, Optional

from swanlab.exceptions import DataStoreError
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest, DeliverSyncStartResponse
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg import executor
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.protocol.core import CoreSyncProtocol

__all__ = ["CoreSyncPython"]


class CoreSyncPython(CoreSyncProtocol):
    def __init__(self):
        super().__init__()
        self._run_ctx: Optional[CoreContext] = None
        self._reader = DataStoreReader()
        self._processor = RecordProcessor()
        self._read_executor = executor.EventLoopExecutor("SwanLab·Sync·Reader")

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
        for record_bytes in self._reader:
            record = Record()
            record.ParseFromString(record_bytes)
            self._processor.process(record)
            # TODO 处理traker记录

    def deliver_sync_finish(self):
        self._read_executor.wait()
        self._reader.close()


class RecordProcessor:
    """
    处理记录的类，根据记录类型调用相应的处理方法
    """

    def __init__(self) -> None:
        # 记录处理函数映射
        self._process_handlers: dict[str, Callable[[Record], None]] = {
            "start": self._process_start,
            "finish": self._process_finish,
            "column": self._process_column,
            "scalar": self._process_scalar,
            "media": self._process_media,
            "config": self._process_config,
            "log": self._process_log,
            "metadata": self._process_metadata,
            "requirements": self._process_requirements,
            "conda": self._process_conda,
            "save": self._process_save,
        }

    def process(self, record: Record):
        """
        处理日志记录
        """
        record_type = record.WhichOneof("record_type")
        handler = self._process_handlers.get(record_type)
        if handler is not None:
            handler(getattr(record, record_type))
        else:
            console.warning(
                f"Skipping unsupported record type during sync: {record_type!r}. "
                "This may indicate that the run file was created by an incompatible SwanLab version."
            )

    def _process_start(self, record: Record): ...

    def _process_finish(self, record: Record): ...

    def _process_column(self, record: Record): ...

    def _process_scalar(self, record: Record): ...

    def _process_media(self, record: Record): ...

    def _process_config(self, record: Record): ...

    def _process_log(self, record: Record): ...

    def _process_metadata(self, record: Record): ...

    def _process_requirements(self, record: Record): ...

    def _process_conda(self, record: Record): ...

    def _process_save(self, record: Record): ...
