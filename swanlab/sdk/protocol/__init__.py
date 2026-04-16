"""
@author: cunyue
@file: __init__.py
@time: 2026/4/12 01:10
@description: SwanLab SDK 协议定义模块
所谓的协议指的是一些可以被用户自定义的接口规范，用户可以通过实现这些协议来扩展 SwanLab SDK 的功能
"""

from .callbacker import Callback
from .core import CoreEnum, CoreProtocol

__all__ = ["Callback", "CoreProtocol", "CoreEnum"]
