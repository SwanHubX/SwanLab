"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:45
@description: SwanLab 运行时上下文，存储Key等存粹的动态信息
此上下文应该只在swanlab.init中使用
"""

from swanlab.sdk.internal.settings import Settings, settings

from .context import RunConfig, RunContext, clear_context, get_context, has_context, set_context, use_temp_context

__all__ = [
    "RunContext",
    "RunConfig",
    "set_context",
    "clear_context",
    "get_context",
    "has_context",
    "use_temp_context",
    "get_current_settings",
]


def get_current_settings() -> Settings:
    """
    获取当前SwanLab配置信息，如果上下文中存在，则使用上下文中的配置信息，否则使用全局配置信息
    """
    if has_context():
        return get_context().config.settings
    return settings
