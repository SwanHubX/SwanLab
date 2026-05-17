"""
@author: cunyue
@file: sync.py
@time: 2026/5/14 14:21
@description: sync 协议实现
"""

from typing import Optional

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.protocol.core import CoreSyncProtocol

__all__ = ["CoreSyncPython"]


class CoreSyncPython(CoreSyncProtocol):
    def __init__(self):
        super().__init__()
        self._run_ctx: Optional[CoreContext] = None

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest):
        pass

    def deliver_sync_finish(self):
        pass
