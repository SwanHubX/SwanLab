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
import os, re


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


def check_key_format(key: str) -> str:
    """检查key字符串格式，必须是0-9a-zA-Z _-和/组成的字符串，并且开头必须是0-9a-zA-Z

    Parameters
    ----------
    key : str
        待检查的字符串
    """
    if not isinstance(key, str):
        raise TypeError(f"key: {key} is not a string")
    # 定义正则表达式
    pattern = re.compile("^[0-9a-zA-Z][0-9a-zA-Z_/-]*$")

    # 检查 key 是否符合规定格式
    if not pattern.match(key):
        raise ValueError(
            f"key: {key} is not a valid string, which must be composed of 0-9a-zA-Z _- and /, and the first character must be 0-9a-zA-Z"
        )
    return key


def check_name_format(name: str, max_len: int = 20) -> str:
    """检查name字符串格式，必须是0-9a-zA-Z _-和/或者中文字符组成的字符串，并且开头必须是0-9a-zA-Z或者中文字符
    最大长度为max_len个字符，一个中文字符算一个字符，如果超出长度，将被截断

    Parameters
    ----------
    name : str
        待检查的字符串

    Returns
    -------
    str
        检查后的字符串
    """
    if not isinstance(name, str):
        raise TypeError(f"name: {name} is not a string")
    # 定义正则表达式
    pattern = re.compile("^[0-9a-zA-Z\u4e00-\u9fa5][0-9a-zA-Z\u4e00-\u9fa5_/-]*$")
    # 检查 name 是否符合规定格式
    if not pattern.match(name):
        raise ValueError(
            f"name: {name} is not a valid string, which must be composed of 0-9a-zA-Z _- and / or Chinese characters, and the first character must be 0-9a-zA-Z or Chinese characters"
        )
    # 检查长度
    if len(name) > max_len:
        name = name[:max_len]
    return name
