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


def check_exp_name_format(name: str, auto_cut: bool = True) -> str:
    """检查实验名格式，必须是0-9a-zA-Z和连字符(_-)，并且不能以连字符(_-)开头或结尾
    最大长度为100个字符，一个中文字符算一个字符

    Parameters
    ----------
    name : str
        待检查的字符串
    auto_cut : bool, optional
        如果超出长度，是否自动截断，默认为True
        如果为False，则超出长度会抛出异常

    Returns
    -------
    str
        检查后的字符串

    Raises
    ------
    TypeError
        name不是字符串，或者name为空字符串
    ValueError
        name不符合规定格式
    IndexError
        name超出长度
    """
    max_len = 100
    if not isinstance(name, str) or name == "":
        raise TypeError(f"name: {name} is not a string")
    # 定义正则表达式
    pattern = re.compile(r"^[0-9a-zA-Z][0-9a-zA-Z_-]*[0-9a-zA-Z]$")
    # 检查 name 是否符合规定格式
    if not pattern.match(name):
        raise ValueError(
            f"name: {name} is not a valid string, which must be composed of 0-9a-zA-Z _- and / or Chinese characters, and the first character must be 0-9a-zA-Z or Chinese characters"
        )
    # 检查长度
    if auto_cut and len(name) > max_len:
        name = name[:max_len]
    elif not auto_cut and len(name) > max_len:
        raise IndexError(f"name: {name} is too long, which must be less than {max_len} characters")
    return name


def check_desc_format(description: str, auto_cut: bool = True):
    """检查实验描述
    不能超过255个字符，可以包含任何字符

    Parameters
    ----------
    description : str
        需要检查和处理的描述信息
    auto_cut : bool
        如果超出长度，是否裁剪并抛弃多余部分

    Returns
    -------
    str
        检查后的字符串，同时会去除字符串头尾的空格

    Raises
    ------
    IndexError
        name超出长度
    """
    max_length = 255
    description = description.strip()

    if len(description) > max_length:
        if auto_cut:
            return description[:max_length]
        else:
            raise IndexError(f"description too long that exceeds {max_length} characters.")
    return description


def check_proj_name_format(name: str, auto_cut: bool = True) -> str:
    """检查项目名格式，必须是0-9a-zA-Z和中文以及连字符(_-)，并且不能以连字符(_-)开头或结尾
    最大长度为100个字符，一个中文字符算一个字符

    Parameters
    ----------
    name : str
        待检查的字符串
    auto_cut : bool, optional
        如果超出长度，是否自动截断，默认为True
        如果为False，则超出长度会抛出异常

    Returns
    -------
    str
        检查后的字符串

    Raises
    ------
    TypeError
        name不是字符串，或者name为空字符串
    ValueError
        name不符合规定格式
    IndexError
        name超出长度
    """
    max_len = 100
    if not isinstance(name, str) or name == "":
        raise TypeError(f"name: {name} is not a string")
    # 定义正则表达式
    pattern = re.compile(r"^[0-9a-zA-Z\u4e00-\u9fa5]+[0-9a-zA-Z\u4e00-\u9fa5_-]*[0-9a-zA-Z\u4e00-\u9fa5]$")
    # 检查 name 是否符合规定格式
    if not pattern.match(name):
        raise ValueError(
            f"name: {name} is not a valid string, which must be composed of 0-9a-zA-Z _- and / or Chinese characters, and the first character must be 0-9a-zA-Z or Chinese characters"
        )
    # 检查长度
    if auto_cut and len(name) > max_len:
        name = name[:max_len]
    elif not auto_cut and len(name) > max_len:
        raise IndexError(f"name: {name} is too long, which must be less than {max_len} characters")
    return name
