from typing import Optional

from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.protocol.core import CoreSyncProtocol

from .context.components.core import CoreBackend, create_core, resolve_core_backend, warn_fallback_once
from .context.components.probe import create_probe

__all__ = ["create_core", "create_probe", "create_core_sync", "resolve_core_sync_backend"]


def resolve_core_sync_backend() -> CoreBackend:
    """解析 sync 实际使用的 core 后端；Go core 二进制不可用时回退内置 Python sync。"""
    return resolve_core_backend("online")


def create_core_sync(backend: Optional[CoreBackend] = None) -> CoreSyncProtocol:
    """core sync 对象工厂；默认走 Go core，缺失二进制时单次告警并回退 CoreSyncPython。"""
    resolved = backend if backend is not None else resolve_core_sync_backend()
    if resolved.impl == "python":
        if resolved.fallback_reason is not None:
            warn_fallback_once(resolved.fallback_reason)
        return CoreSyncPython()

    from swanlab.sdk.internal.core_go import CoreSyncGo

    return CoreSyncGo(binary=resolved.binary)
