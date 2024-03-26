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
import ujson
import yaml


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


def check_string(target: str):
    """
    不能全空格，也不能为空字符串
    """
    # 利用正则表达式匹配非空格字符
    if re.match(r"^\s*$", target):
        return False
    # 利用正则表达式匹配非空字符串
    if re.match(r"^\s*$", target) or target == "":
        return False
    return True


def check_tag_format(key: str, auto_cut=True) -> str:
    """检查tag字符串格式，必须是0-9a-zA-Z _-和/组成的字符串(包含空格)，并且开头必须是0-9a-zA-Z
    最大长度为255字符

    Parameters
    ----------
    key : str
        待检查的字符串
    """
    max_len = 255
    if not isinstance(key, str):
        raise TypeError(f"tag: {key} is not a string")
    if not check_string(key):
        raise ValueError(f"tag: {key} is an empty string")
    # 检查长度
    if auto_cut and len(key) > max_len:
        key = key[:max_len]
    elif not auto_cut and len(key) > max_len:
        raise IndexError(f"tag: {key} is too long, which must be less than {max_len} characters")
    return key.strip()


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
    if not check_string(name):
        raise ValueError(f"name: {name} is an empty string")
    # 检查长度
    if auto_cut and len(name) > max_len:
        name = name[:max_len]
    elif not auto_cut and len(name) > max_len:
        raise IndexError(f"name: {name} is too long, which must be less than {max_len} characters")
    return name.strip()


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
    if not check_string(name):
        raise ValueError(f"name: {name} is an empty string")
    # 检查长度
    if auto_cut and len(name) > max_len:
        name = name[:max_len]
    elif not auto_cut and len(name) > max_len:
        raise IndexError(f"name: {name} is too long, which must be less than {max_len} characters")
    return name.strip()


def check_load_json_yaml(file_path: str, paramter_name: str = "init_path"):
    # 不是字符串
    if not isinstance(file_path, str):
        raise TypeError("{} must be a string, but got {}".format(paramter_name, type(file_path)))
    # 检查file_path的后缀是否是json/yaml，否则报错
    path_suffix = file_path.split(".")[-1]
    if not file_path.endswith((".json", ".yaml", ".yml")):
        raise ValueError(
            "{} must be a json or yaml file ('.json', '.yaml', '.yml'), but got {}, please check if the content of config_file is correct.".format(
                paramter_name, path_suffix
            )
        )
    # 转换为绝对路径
    file_path = os.path.abspath(file_path)
    # 读取配置文件
    # 如果文件不存在或者不是文件
    if (not os.path.exists(file_path)) or (not os.path.isfile(file_path)):
        raise FileNotFoundError("{} not found, please check if the file exists.".format(paramter_name))
    # 为空
    if os.path.getsize(file_path) == 0:
        raise ValueError("{} is empty, please check if the content of config_file is correct.".format(paramter_name))
    # 无权限读取
    if not os.access(file_path, os.R_OK):
        raise PermissionError(
            "No permission to read {}, please check if you have the permission.".format(paramter_name)
        )
    load = ujson.load if path_suffix == "json" else yaml.safe_load
    with open(file_path, "r") as f:
        # 读取配置文件的内容
        file_data = load(f)
        # 如果读取的内容不是字典类型，则报错
        if not isinstance(file_data, dict):
            raise TypeError("The configuration file must be a dictionary, but got {}".format(type(file_data)))
    return file_data


# ---------------------------------- 一些格式检查的工具函数 ----------------------------------


def is_ipv4(string: str) -> bool:
    """判断字符串是否是一个ipv4地址

    Parameters
    ----------
    string : str
        待检查的字符串

    Returns
    -------
    bool
        如果是ipv4地址，返回True，否则返回False
    """
    pattern = re.compile(r"^((2[0-4]\d|25[0-5]|[01]?\d\d?)\.){3}(2[0-4]\d|25[0-5]|[01]?\d\d?)$")
    return isinstance(string, str) and pattern.match(string)


def is_port(string: str) -> bool:
    """判断字符串是否是一个端口号

    Parameters
    ----------
    string : str
        待检查的字符串

    Returns
    -------
    bool
        如果是端口号，返回True，否则返回False
    """
    if not is_int(string):
        return False
    port = int(string)
    return 0 <= port <= 65535


def is_int(string: str) -> bool:
    """判断字符串是否可以转换为整数

    Parameters
    ----------
    string : str
        待检查的字符串

    Returns
    -------
    bool
        如果可以转换为整数，返回True，否则返回False
    """
    try:
        int(string)
        return True
    except ValueError:
        return False
