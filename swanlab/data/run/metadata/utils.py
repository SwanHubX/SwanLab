"""
@author: cunyue
@file: utils.py
@time: 2025/9/9 12:38
@description: 一些元信息采集工具函数
"""

import os

import wrapt


def check_env(key: str, default_return=None):
    """
    函数执行收尾，类似一个开关，判断环境变量中的对应 key 是否存在并符合预期值，如果是则不执行函数，直接返回默认值
    :param key: 环境变量中的对应 key
    :param default_return: 如果环境变量中的对应 key 存在并符合预期值，则返回的默认值
    :return: 装饰器
    """

    @wrapt.decorator
    def wrapper(wrapped, _, args, kwargs):
        disable = os.getenv(key, "").lower()
        if disable in ["true", "1", "yes", "on"]:
            return default_return
        return wrapped(*args, **kwargs)

    return wrapper
