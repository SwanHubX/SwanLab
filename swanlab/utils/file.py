#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-02 13:23:42
@File: swanlab\utils\file.py
@IDE: vscode
@Description:
    文件操作
"""
import portalocker
from functools import wraps
from io import TextIOWrapper
import os


# 锁定文件，防止多进程写入同一个文件
# 这是一个装饰器，用于锁定文件，
def lock_file(file_path: str, mode: str = "r+"):
    """锁定文件，防止多进程写入同一个文件
    装饰器将向函数传递一个file参数，这个参数是一个文件对象，可以直接写入
    不需要手动打开和关闭文件，否则导致锁定失效

    Parameters
    ----------
    file_path : str
        文件路径
    mode : str, optional
        文件读写模式, by default "r+"
    """

    def decorator(func):
        # 保证函数签名不变
        @wraps(func)
        def wrapper(*args, **kwargs):
            f = open(file_path, mode=mode, encoding="utf-8")
            portalocker.lock(f, portalocker.LOCK_EX)
            try:
                res = func(file=f, *args, **kwargs)
            except Exception as e:
                raise e
            finally:
                portalocker.unlock(f)
                f.close()
            return res

        return wrapper

    return decorator


def get_a_lock(file_path: str, mode: str = "r+", encoding="utf-8") -> TextIOWrapper:
    """获取一个文件锁,
    返回文件对象，你需要手动关闭文件
    """
    f = open(file_path, mode=mode, encoding=encoding)
    portalocker.lock(f, portalocker.LOCK_EX)
    return f
