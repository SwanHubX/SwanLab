"""
@author: cunyue
@file: degraded.py
@time: 2026/7/30 16:41
@description: Settings 初始化失败时使用的降级代理。
"""

from typing import NoReturn


class DegradedSettings:
    """全局降级单例代理：配置访问会抛出原始初始化错误。"""

    def __init__(self, error: BaseException) -> None:
        object.__setattr__(self, "_degraded_error", error)

    def raise_for_access(self, name: str) -> NoReturn:
        err = object.__getattribute__(self, "_degraded_error")
        raise RuntimeError(
            f"SwanLab settings failed to initialize; cannot access `{name}`. "
            "Fix the configuration and restart the process."
        ) from err

    def __getattr__(self, name: str) -> NoReturn:
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        self.raise_for_access(name)

    def __repr__(self) -> str:
        err = object.__getattribute__(self, "_degraded_error")
        return f"<DegradedSettings: {type(err).__name__}: {err}>"
