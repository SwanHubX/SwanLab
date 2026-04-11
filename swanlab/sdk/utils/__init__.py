"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:10
@description: 可供外部使用的 SwanLab SDK 工具包
理论上本模块的内容都可以被用户调用，我们在此处导出一些常用函数，并将写入API文档中
"""

from swanlab.sdk.internal.pkg.helper import get_swanlab_latest_version
from swanlab.sdk.internal.pkg.helper.env import DEBUG

__all__ = ["DEBUG", "get_swanlab_latest_version"]
