"""
@author: cunyue
@file: __init__.py
@time: 2026/7/30 16:41
@description: 设置门控：控制构造 Settings 时是否加载外部配置源。

门控信号通过线程本地存储承载，避免类级可变状态在多线程并发构造时跨线程泄漏。
仅用于 import 期间的降级构造——环境变量非法时，构造一个不读任何外部配置的纯默认实例，
避免在 import swanlab 时直接抛错；真正的校验错误会在后续 swanlab.init() 新建 Settings() 时照常抛出。


由于各子模块（core / experiment / integration 等）的向下兼容环境变量绕过了pydantic的加载机制，
因此读取统一通过 env_or 接入门控，保证降级实例为不含任何外部配置的纯默认值。

作为工具模块，为了区分于配置模块，改用文件夹申明模块
"""

import os
import threading
from typing import Optional, overload


class _SettingsLocal(threading.local):
    """线程本地门控状态：load_external 为 True 时正常加载外部配置源。"""

    load_external: bool = True


settings_local = _SettingsLocal()


def load_external_enabled() -> bool:
    """读取当前线程的门控状态，默认 True（加载外部配置源）。"""
    return settings_local.load_external


@overload
def env_or(key: str, default: str) -> str: ...


@overload
def env_or(key: str, default: None) -> Optional[str]: ...


def env_or(key: str, default: Optional[str]) -> Optional[str]:
    """读取环境变量；降级构造时（门控关闭）直接返回 default，不触碰 env。

    用于向下兼容的旧版环境变量读取：保证降级实例为不含任何外部配置的纯默认值。
    """
    if not load_external_enabled():
        return default
    return os.environ.get(key, default)
