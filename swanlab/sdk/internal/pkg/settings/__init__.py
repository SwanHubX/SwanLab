"""
@author: cunyue
@file: __init__.py
@time: 2026/3/8 15:41
@description: SwanLab SDK 配置信息工具
"""

from swanlab.sdk.internal.context import get_context, has_context
from swanlab.sdk.internal.settings import Settings, settings


def get_current_settings() -> Settings:
    """
    获取当前SwanLab配置信息，如果上下文中存在，则使用上下文中的配置信息，否则使用全局配置信息
    """
    if has_context():
        return get_context().config.settings
    return settings
