from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.protocol.core import CoreSyncProtocol

from .context.components.core import create_core
from .context.components.probe import create_probe

__all__ = ["create_core", "create_probe", "create_core_sync"]


def create_core_sync() -> CoreSyncProtocol:
    if helper.get_core_impl() == "python":
        return CoreSyncPython()
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError("The SwanLab Go core sync runtime is not available yet.")
