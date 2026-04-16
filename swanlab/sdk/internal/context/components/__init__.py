"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 21:12
@description: 上下文组件，存储swanlab run 运行生命周期所需要的外置组件，他们包括：

1. core: 核心组件，负责处理record的持久化和后端交互
2. system: 系统组件，负责系统环境信息采集、监控指标上报
3. callback: 回调组件，负责处理用户自定义的回调函数
4. record_builder: 记录构建器组件，将事件转换为record

PS: 考虑到未来多语言SDK的支持，在设计上core和system在未来会被当作独立二进制组件实现

TODO system 组件暂时挂在 run 对象上，后续迁移到上下文中
"""

from .callbacker import CallbackManager, callbacker, create_callback_manager
from .core import create_core

__all__ = ["CallbackManager", "callbacker", "create_callback_manager", "create_core"]
