#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
    除了utils和error模块，其他模块都可以使用这个模块
"""
import enum
import os
from typing import List

import swankit.env as E
from swankit.env import SwanLabSharedEnv


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
    MODE = SwanLabSharedEnv.SWANLAB_MODE.value
    """
    swanlab的解析模式，涉及操作员注册的回调，目前有三种：local、cloud、disabled，默认为cloud
    大小写敏感
    """
    SWANBOARD_PROT = "SWANLAB_BOARD_PORT"
    """
    cli swanboard 服务端口
    """
    SWANBOARD_HOST = "SWANLAB_BOARD_HOST"
    """
    cli swanboard 服务地址
    """
    WEB_HOST = "SWANLAB_WEB_HOST"
    """
    swanlab云端环境的web地址
    """
    API_HOST = "SWANLAB_API_HOST"
    """
    swanlab云端环境的api地址
    """
    API_KEY = "SWANLAB_API_KEY"
    """
    云端api key，登录时会首先查找此环境变量，如果不存在，判断用户是否已登录，未登录则进入登录流程

    * 如果login接口传入字符串，此环境变量无效，此时相当于绕过 get_key 接口
    * 如果用户已登录，此环境变量的优先级高于本地存储登录信息
    """
    RUNTIME = "SWANLAB_RUNTIME"
    """
    swanlab的运行时环境，"user" "develop" "test" "test-no-cloud" "task"
    """
    WEBHOOK = "SWANLAB_WEBHOOK"
    """
    webhook地址。swanlab初始化完毕时，如果此环境变量存在，会调用此地址，发送消息。
    """

    @classmethod
    def set_default(cls):
        """
        设置默认的环境变量值
        """
        envs = {
            cls.WEB_HOST.value: "https://swanlab.cn",
            cls.API_HOST.value: "https://api.swanlab.cn/api",
            cls.RUNTIME.value: "user",
        }
        for k, v in envs.items():
            os.environ.setdefault(k, v)

    @classmethod
    def check(cls):
        """
        检查环境变量的值是否为预期值中的一个
        :raises ValueError: 如果环境变量的值不在预期值中
        """
        envs = {
            cls.MODE.value: ["local", "cloud", "disabled"],
            cls.RUNTIME.value: ["user", "develop", "test", "test-no-cloud", "task"],
        }
        for k, vs in envs.items():
            if k in os.environ and os.environ[k] not in vs:
                raise ValueError(f"Unknown value for {k}: {os.environ[k]}")

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
