"""
@author: cunyue
@file: __init__.py.py
@time: 2025/6/16 13:21
@description: swanlab 核心业务代码 - python 版，将包含：
1. 上传线程逻辑
2. http 客户端代码
"""

# FIXME 存在循环引用，我们需要更优雅的代码结构
# from . import auth
# from . import uploader
from .client import *
from .utils import timer
