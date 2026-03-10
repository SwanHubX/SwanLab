"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:10
@description: 可供外部使用的 SwanLab SDK 工具包
理论上本模块的内容都可以被用户调用，我们在此处导出一些常用函数，并将写入API文档中
"""

from .callbacker import SwanLabCallback
from .experiment import generate_color, generate_id, generate_name
from .helper.env import DEBUG
from .version import get_swanlab_latest_version

__all__ = ["generate_color", "generate_id", "generate_name", "DEBUG", "get_swanlab_latest_version", "SwanLabCallback"]
