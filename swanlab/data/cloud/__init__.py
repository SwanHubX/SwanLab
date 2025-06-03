#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 16:42
@File: __init__.py
@IDE: pycharm
@Description:
    云端日志资源上传部分
    NOTE: 这部分设计已经过时，后续考虑优化
"""
from .start_thread import ThreadPool
from .task_types import UploadType

__all__ = ["UploadType", "ThreadPool"]
