"""
@author: cunyue
@file: rich.py
@time: 2026/3/8 12:12
@description: 封装rich包的一些常用功能
"""

from functools import wraps

from rich.status import Status

__all__ = ["with_loading_animation"]


def with_loading_animation(message: str = "Initializing SwanLab...", spinner_name: str = "dots"):
    """
    一个使用 rich 包装的加载动画装饰器
    :param message: 动画旁边显示的文字
    :param spinner_name: 动画样式（如 'dots', 'bouncingBar', 'point' 等）
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # 使用 rich 的 status 作为上下文管理器包裹函数的执行
            with Status(message, spinner=spinner_name):
                return func(*args, **kwargs)

        return wrapper

    return decorator
