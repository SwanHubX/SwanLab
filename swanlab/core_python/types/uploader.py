"""
@author: cunyue
@file: uploader.py
@time: 2026/3/3 17:38
@description: 上传器类型定义
"""

from typing import Callable


UploadCallback = Callable[[int], None]
"""
上传进度回调函数，参数为本次已上传的数量
"""
