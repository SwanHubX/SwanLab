#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
    除了utils和error模块，其他模块都可以使用这个模块
"""
import os
from typing import List
import swankit.env as E
from swankit.env import SwanLabSharedEnv
import enum


# ---------------------------------- 环境变量枚举类 ----------------------------------


class SwanLabEnv(enum.Enum):
    """
    swanlab环境变量枚举类，包含swankit的共享环境变量
    """
    SWANLAB_FOLDER = SwanLabSharedEnv.SWANLAB_FOLDER.value
    """
    swanlab全局文件夹保存的路径，默认为用户主目录下的.swanlab文件夹
    """
    SWANLOG_FOLDER = SwanLabSharedEnv.SWANLOG_FOLDER.value
    """
    swanlab解析日志文件保存的路径，默认为当前运行目录的swanlog文件夹
    """
    SWANLAB_MODE = SwanLabSharedEnv.SWANLAB_MODE.value
    """
    swanlab的解析模式，涉及操作员注册的回调，目前有三种：local、cloud、disabled，默认为cloud
    大小写不敏感
    """
    SWANLAB_PORT = "SWANLAB_SERVER_PORT"
    """cli 服务端口"""
    SWANLAB_HOST = "SWANLAB_SERVER_HOST"
    """cli 服务地址"""
    SWANLAB_PACKAGE = "SWANLAB_PACKAGE_PATH"
    """
    swanlab的包路径，即package.json文件路径，可以设置为相对路径，但最终会转换为绝对路径
    """

    @classmethod
    def list(cls) -> List[str]:
        """
        获取所有的枚举值
        :return: 所有的枚举值
        """
        return [item.value for item in cls]


# ---------------------------------- API ----------------------------------

is_windows = E.is_windows

get_mode = E.get_mode

get_swanlog_dir = E.get_swanlog_dir

get_save_dir = E.get_save_dir


def in_jupyter() -> bool:
    """
    用于检测是否在 notebook jupyter 中运行
    :return: bool 是否在 notebook 中运行
    """
    try:
        # notebook 中会有 __IPYTHON__，而正常环境没有定义，所以 try
        # 'type: ignore': 可以让 pylance 忽略对变量定义的检查
        _ = __IPYTHON__  # type: ignore
        return True
    except NameError:
        return False


def get_package_path() -> str:
    """
    获取swanlab的包路径，即package.json文件路径
    :raise FileNotFoundError: 文件不存在时抛出异常
    :raise IsADirectoryError: 文件是一个文件夹时抛出异常
    :return: swanlab的包路径，是一个绝对路径
    """
    path = os.getenv(SwanLabEnv.SWANLAB_PACKAGE.value) or os.path.join(os.path.dirname(__file__), "package.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"package.json not found in {path}")
    if os.path.isdir(path):
        raise IsADirectoryError(f"{path} is a directory")
    return os.path.abspath(path)

############################################
# port和host不会在这里设置，而是在cli模块中设置
############################################
