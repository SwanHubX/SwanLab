#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 16:39:40
@File: swanlab/server/module/utils.py
@IDE: vscode
@Description:
    api的一些工具函数
"""
from functools import wraps
from .resp import UNEXPECTED_ERROR_500


def catch_error():
    """捕获异常的装饰器"""

    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                return await func(*args, **kwargs)
            except Exception as e:
                return UNEXPECTED_ERROR_500(str(e))

        return wrapper

    return decorator
