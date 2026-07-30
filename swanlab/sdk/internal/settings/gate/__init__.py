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


# --------------------------------------------------------------------------------------------------
# 降级状态追踪
#
# 全局单例 settings 在 import 期降级构造时，记住导致降级的原始异常。
# 消费降级单例中 host / api_key 等不可靠值的入口（Api、login）通过 raise_if_degraded() re-raise，
# 避免因降级默认值导致请求被路由到公有云或抛出误导性错误。
#
# 仅在 import 期写入一次，进程内只读，无需线程本地存储。
# --------------------------------------------------------------------------------------------------
_degraded_error: Optional[BaseException] = None


def set_degraded(error: BaseException) -> None:
    """记录降级构造时的原始异常。仅在全局单例降级时调用一次。"""
    global _degraded_error
    _degraded_error = error


def raise_if_degraded() -> None:
    """全局单例处于降级状态时 re-raise 原始配置错误。

    供消费降级单例中 host / api_key 等不可靠值的入口调用，
    让真正的配置错误在使用时暴露，而非静默路由到公有云或抛出误导性错误。
    """
    if _degraded_error is not None:
        raise _degraded_error


class DegradedSettings:
    """全局降级单例代理：任意属性访问 re-raise 原始配置错误。

    替代在各消费入口手动埋点：降级单例本身即为"不可信"标记，
    任何对 host / api_key 等字段的读取都会立即暴露真实配置错误，
    避免请求被路由到公有云或抛出误导性错误。

    显式参数场景天然正确——调用方传入了完整的 api_key + host 时不会触碰单例，因此不触发 __getattr__。
    """

    def __init__(self, error: BaseException) -> None:
        object.__setattr__(self, "_degraded_error", error)

    def __getattr__(self, name: str):
        raise object.__getattribute__(self, "_degraded_error")

    def __repr__(self) -> str:
        err = object.__getattribute__(self, "_degraded_error")
        return f"<DegradedSettings: {type(err).__name__}: {err}>"
