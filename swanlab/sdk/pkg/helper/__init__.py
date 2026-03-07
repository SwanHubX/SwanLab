"""
@author: cunyue
@file: __init__.py.py
@time: 2026/3/7 14:18
@description: SwanLab SDK 辅助函数
"""

from functools import wraps


def catch_and_return_none(*exceptions):
    """
    Catch the specified exception and return None.
    If no argument is passed, it defaults to catching Exception (though this is not recommended).
    """
    # 如果没传具体的异常类，就用最宽泛的 Exception 兜底
    catch_types = exceptions if exceptions else (Exception,)

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except catch_types:
                return None

        return wrapper

    return decorator
