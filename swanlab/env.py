#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
    除了utils和error模块，其他模块都可以使用这个模块
"""
import datetime
import enum
import io
import netrc
import os
import re
import sys
from enum import Enum
from pathlib import Path
from platform import platform
from typing import List, Union
from urllib.parse import urlparse


class SwanLabMode(Enum):
    """
    swanlab的解析模式，枚举类
    """

    DISABLED = "disabled"
    CLOUD = "cloud"
    OFFLINE = "offline"
    LOCAL = "local"

    @classmethod
    def list(cls) -> List[str]:
        """
        获取所有的枚举值
        :return: 所有的枚举值
        """
        return [item.value for item in cls]


class SwanLabEnv(enum.Enum):
    """
    swanlab环境变量枚举类，包含swankit的共享环境变量
    """

    LOG_LEVEL = "SWANLAB_LOG_LEVEL"
    """
    swanlab日志的输出级别，默认为info，可选值有debug、info、warning、error、critical
    """
    SWANLAB_FOLDER = "SWANLAB_SAVE_DIR"
    """
    swanlab全局文件夹保存的路径，默认为用户主目录下的.swanlab文件夹
    """
    SWANLOG_FOLDER = "SWANLAB_LOG_DIR"
    """
    swanlab解析日志文件保存的路径，默认为当前运行目录的swanlog文件夹
    """
    MODE = "SWANLAB_MODE"
    """
    swanlab的解析模式，涉及操作员注册的回调，目前有四种：local、cloud、disabled、offline，默认为cloud
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
    WORKSPACE = "SWANLAB_WORKSPACE"
    """
    swanlab的工作空间，默认为当前登录用户
    """
    PROJ_NAME = "SWANLAB_PROJ_NAME"
    """
    swanlab的项目名称
    """
    EXP_NAME = "SWANLAB_EXP_NAME"
    """
    swanlab的实验名称
    """
    RUN_ID = "SWANLAB_RUN_ID"
    """
    swanlab 实验运行 id，resume 时使用
    """
    RESUME = "SWANLAB_RESUME"
    """
    swanlab 实验是否为恢复运行，resume 时使用，设置为 true 则表示恢复运行
    """
    RUNTIME = "SWANLAB_RUNTIME"
    """
    swanlab的运行时环境，"user" "develop" "test" "test-no-cloud"，目前仅用于开发测试
    """
    WEBHOOK = "SWANLAB_WEBHOOK"
    """
    webhook地址。swanlab初始化完毕时，如果此环境变量存在，会调用此地址，发送消息。
    """
    WEBHOOK_VALUE = "SWANLAB_WEBHOOK_VALUE"
    """
    webhook发送的自定义内容，以字符串读取此环境变量值后发送
    """
    DESCRIPTION = "SWANLAB_DESCRIPTION"
    """
    实验描述，用于为实验提供更详细的介绍或标注
    """
    JOB = "SWANLAB_JOB_TYPE"
    """
    实验任务类型，用于标注当前实验的任务类型，例如分类、回归等
    """
    GROUP = "SWANLAB_GROUP"
    """
    实验组，用于将实验划分到不同的组别，便于管理和区分
    """
    TAGS = "SWANLAB_TAGS"
    """
    实验标签，用于标注当前实验，多个标签用逗号分隔
    """
    DISABLE_GIT = "SWANLAB_DISABLE_GIT"
    """
    禁用Git功能，设置为true时不会采集Git信息
    """

    @staticmethod
    def is_hostname(value: str) -> bool:
        """
        判断是否为合法的主机名（支持 http/https 协议和端口号）
        :param value: 待判断的字符串
        :return: 是否为合法的主机名
        """
        # 解析 URL 以提取主机部分
        parsed_url = urlparse(value)
        if parsed_url.scheme in ["http", "https"]:
            value = parsed_url.hostname  # 只取域名部分

        # 处理 IP 地址
        if re.fullmatch(r'((25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})(\.|$)){4}', value):
            return True

        # 处理域名（去掉端口）
        if re.fullmatch(r'^(?!-)([a-zA-Z0-9-]{1,63}\.)+[a-zA-Z]{2,63}$', value):
            return True

        return False

    @classmethod
    def set_default(cls):
        """
        设置默认的环境变量值，如果netrc文件存在，使用netrc中第一个machine的信息
        """
        # 从netrc文件中读取host信息
        path = os.path.join(get_save_dir(), ".netrc")
        web_host = "https://swanlab.cn"
        api_host = "https://api.swanlab.cn/api"
        if os.path.exists(path):
            nrc = netrc.netrc(path)
            for host, info in nrc.hosts.items():
                base_host = remove_host_suffix(host, "/api")
                web_host = info[0] if cls.is_hostname(info[0]) else base_host
                api_host = base_host + "/api"
                break

        envs = {
            cls.WEB_HOST.value: web_host,
            cls.API_HOST.value: api_host,
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
            cls.MODE.value: ["local", "cloud", "disabled", "offline"],
            cls.RUNTIME.value: ["user", "develop", "test", "test-no-cloud"],
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


def utc_time() -> datetime.datetime:
    """获取当前时间(UTC时区)"""
    return datetime.datetime.now(datetime.timezone.utc)


def create_time() -> str:
    """获取当前时间的iso字符串(UTC时区)"""
    return utc_time().isoformat()


def is_windows() -> bool:
    """判断当前操作系统是否是windows还是类unix系统，主要是路径分隔上的差别
    此外的系统会报错为 UnKnownSystemError
    :raise OSError: 未知系统错误，此时swanlab运行在未知系统上，这个系统不是windows或者类unix系统
    :return: True表示是windows系统，False表示是类unix系统
    """
    if sys.platform.startswith("win"):
        return True
    elif sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
        return False
    raise OSError("Unknown system, not windows or unix-like system")


def is_macos() -> bool:
    """判断当前操作系统是否是macos
    :return: True表示是macos系统，False表示不是macos系统
    """
    return "mac" in platform().lower()


def get_mode() -> str:
    """
    获取当前的swanlab解析模式，如果没有设置，默认为cloud
    :raise ValueError: 未知的swanlab模式
    :return: swanlab的解析模式
    """
    mode = os.getenv(SwanLabEnv.MODE.value)
    if mode is None:
        mode = SwanLabMode.CLOUD.value
    mode = mode.lower()
    if mode not in SwanLabMode.list():
        raise ValueError(f"Unknown swanlab mode: {mode}")
    return mode


def get_save_dir() -> str:
    """
    获取存放swanlab全局文件的文件夹路径，如果不存在就创建
    此函数对应为SWANLAB_SAVE_FOLDER全局变量，如果没有设置，默认为用户主目录下的.swanlab文件夹
    执行此函数时，如果文件夹不存在，自动创建，但是出于安全考虑，不会自动创建父文件夹
    :raises
        :raise FileNotFoundError: folder的父目录不存在
        :raise NotADirectoryError: folder不是一个文件夹
    :return: swanlab全局文件夹保存的路径，返回处理后的绝对路径
    """
    folder = os.getenv(SwanLabEnv.SWANLAB_FOLDER.value)
    if folder is None:
        folder = os.path.join(os.path.expanduser("~"), ".swanlab")
    folder = os.path.abspath(folder)
    if not os.path.exists(os.path.dirname(folder)):
        raise FileNotFoundError(f"{os.path.dirname(folder)} not found")
    if not os.path.exists(folder):
        # 只创建当前文件夹，不创建父文件夹，之所以还要捕捉 FileExistsError，是因为在多线程或多进程环境下，可能会有多个线程或进程同时创建同一个文件夹
        # 比如：https://github.com/SwanHubX/SwanLab/issues/1033
        try:
            os.mkdir(folder)
        except FileExistsError:
            pass
    if not os.path.isdir(folder):
        raise NotADirectoryError(f"{folder} is not a directory")
    return folder


def get_swanlog_dir() -> str:
    """
    获取存放swanlog日志文件的文件夹路径
    此函数对应为SWANLAB_LOG_FOLDER全局变量，如果没有设置，默认为当前运行目录下的swanlog文件夹
    需要注意，此函数并不会保证文件夹的存在，但是会检查父文件夹是否存在以及folder是否是一个文件夹
    :raises
        :raise FileNotFoundError: folder的父目录不存在
        :raise NotADirectoryError: folder不是一个文件夹
    :return: swanlog日志文件保存的路径，返回处理后的绝对路径
    """
    folder = os.getenv(SwanLabEnv.SWANLOG_FOLDER.value)
    if folder is None:
        folder = os.path.join(os.getcwd(), "swanlog")
    folder = os.path.abspath(folder)
    if not os.path.exists(os.path.dirname(folder)):
        raise FileNotFoundError(f"{os.path.dirname(folder)} not found")
    if not os.path.exists(folder):
        return folder
    if not os.path.isdir(folder):
        raise NotADirectoryError(f"{folder} is not a directory")
    return folder


def create_swanlog_dir(logdir: Union[Path, str] = None):
    """
    获取swanlog文件夹，如果文件夹不存在则创建
    :param logdir: swanlog文件夹路径，默认为当前工作目录下的swanlog文件夹
    """
    if logdir is None:
        logdir = get_swanlog_dir()
    try:
        os.makedirs(logdir, exist_ok=True)
        if not os.access(logdir, os.W_OK):
            raise IOError(f"no write permission for path: {logdir}")
    except Exception as error:
        raise IOError(f"Failed to create or access logdir: {logdir}, error: {error}")
    # 如果logdir是空的，创建.gitignore文件，写入*
    if not os.listdir(logdir):
        with open(os.path.join(logdir, ".gitignore"), "w", encoding="utf-8") as f:
            f.write("*")
    return logdir


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


def is_interactive():
    """
    是否为可交互式环境（输入连接tty设备）
    特殊的环境：jupyter notebook
    """
    try:
        fd = sys.stdin.fileno()
        return os.isatty(fd) or in_jupyter()
    # 当使用capsys、capfd或monkeypatch等fixture来捕获或修改标准输入输出时，会抛出io.UnsupportedOperation
    # 多为测试情况，可交互
    except io.UnsupportedOperation:
        return True


def remove_host_suffix(host: str, suffix: str) -> str:
    """
    移除host的后缀
    :param host: 待处理的host
    :param suffix: 要移除的后缀
    :return: 处理后的host
    """
    host = host.rstrip()
    if len(suffix) == 0:
        return host
    if host.endswith(suffix):
        return host[: -len(suffix)]
    return host
