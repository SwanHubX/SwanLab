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
from typing import MutableMapping, Optional
from .utils.file import is_port, is_ipv4
from .error import UnKnownSystemError
import sys

Env = Optional[MutableMapping]

_env = dict()
"""运行时环境变量参数存储，实际上就是一个字典"""

# '描述' = "key"
# ---------------------------------- 基础环境变量 ----------------------------------
ROOT = "SWANLAB_LOG_DIR"
"""命令执行目录SWANLAB_LOG_DIR，日志文件存放在这个目录下，如果自动生成，则最后的目录名为swanlog，第一次调用时如果路径不存在，会自动创建路径"""

PORT = "SWANLAB_SERVER_PORT"
"""服务端口SWANLAB_SERVER_PORT，服务端口"""

HOST = "SWANLAB_SERVER_HOST"
"""服务端口SWANLAB_SERVER_PORT，服务地址"""

DEV = "SWANLAB_DEV"


def get_swanlog_dir(env: Optional[Env] = None) -> Optional[str]:
    """获取swanlog路径

    Returns
    -------
    Optional[str]
        swanlog目录路径
    """
    if _env.get(ROOT) is not None:
        return _env.get(ROOT)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    # 默认为当前目录下的swanlog目录
    default: Optional[str] = os.path.join(os.getcwd(), "swanlog")
    path = env.get(ROOT, default=default)
    # 必须是一个绝对路径
    if not os.path.isabs(path):
        raise ValueError('SWANLAB_LOG_DIR must be an absolute path, now is "{path}"'.format(path=path))
    # 路径必须存在
    if not os.path.exists(path):
        if path == default:
            raise ValueError(
                'The log file was not found in the default path "{path}". Please use the "swanlab watch -l <LOG '
                'PATH>" command to specify the location of the log path."'.format(path=path)
            )
        else:
            raise ValueError('SWANLAB_LOG_DIR must be an existing path, now is "{path}"'.format(path=path))
    # 路径必须是一个目录
    if not os.path.isdir(path):
        raise ValueError('SWANLAB_LOG_DIR must be a directory, now is "{path}"'.format(path=path))
    _env[ROOT] = path
    return path


def get_server_port(env: Optional[Env] = None) -> Optional[int]:
    """获取服务端口

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map,可以是任意实现了MutableMapping的对象, 默认将使用os.environ

    Returns
    -------
    Optional[int]
        服务端口
    """
    # 第一次调用时，从环境变量中提取，之后就不再提取，而是从缓存中提取
    if _env.get(PORT) is not None:
        return _env.get(PORT)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    default: Optional[int] = 5092
    port = env.get(PORT, default=default)
    # 必须可以转换为整数，且在0-65535之间
    if not is_port(port):
        raise ValueError('SWANLAB_SERVER_PORT must be a port, now is "{port}"'.format(port=port))
    _env[PORT] = int(port)
    return _env.get(PORT)


def get_server_host(env: Optional[Env] = None) -> Optional[str]:
    """获取服务端口

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map,可以是任意实现了MutableMapping的对象, 默认将使用os.environ

    Returns
    -------
    Optional[int]
        服务端口
    """
    default: Optional[str] = "127.0.0.1"
    # 第一次调用时，从环境变量中提取，之后就不再提取，而是从缓存中提取
    if _env.get(HOST) is not None:
        return _env.get(HOST)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    _env[HOST] = env.get(HOST, default=default)
    # 必须是一个ipv4地址
    if not is_ipv4(_env.get(HOST)):
        raise ValueError('SWANLAB_SERVER_HOST must be an ipv4 address, now is "{host}"'.format(host=_env.get(HOST)))
    return _env.get(HOST)


def is_dev(env: Optional[Env] = None) -> bool:
    """判断是否是开发模式

    Returns
    -------
    bool
        是否是开发模式
    """
    if _env.get(DEV) is not None:
        return _env.get(DEV) == "TRUE"
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    _env[DEV] = env.get(DEV, default=False)
    return _env.get(DEV) == "TRUE"


# ---------------------------------- 初始化基础环境变量 ----------------------------------

# 所有的初始化函数
function_list = [
    get_swanlog_dir,
    get_server_port,
    get_server_host,
    is_dev
]


# 定义初始化函数
def init_env(env: Optional[Env] = None):
    """初始化环境变量

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map,可以是任意实现了MutableMapping的对象, 默认将使用os.environ
    """
    for func in function_list:
        func(env)


# ---------------------------------- 计算变量 ----------------------------------
DATABASE_PATH = "SWANLAB_DB_PATH"
"""日志目录SWANLAB_LOG_DIR，日志文件存放在这个目录下"""


# ---------------------------------- 定义变量访问方法 ----------------------------------


def get_db_path() -> Optional[str]:
    """获取数据库路径，这是一个计算变量，
    通过`get_swanlog_dir()`返回值得到

    Returns
    -------
    Optional[str]
        数据库文件路径
    """
    if _env.get(DATABASE_PATH) is not None:
        return _env.get(DATABASE_PATH)
    # 否则从环境变量中提取
    _env[DATABASE_PATH] = os.path.join(get_swanlog_dir(), "runs.swanlab")
    return _env.get(DATABASE_PATH)


def is_windows() -> bool:
    """判断当前操作系统是否是windows还是类unix系统
    此外的系统会报错为 UnKnownSystemError

    Returns
    -------
    bool
        是否是windows
    """
    if sys.platform.startswith("win"):
        return True
    elif sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
        return False
    raise UnKnownSystemError("Unknown system, not windows or unix-like system")


def get_user_home() -> str:
    """获取用户家目录，需要分为windows和类unix系统

    Returns
    -------
    str
        用户家目录
    """
    if is_windows():
        return os.environ.get("USERPROFILE")
    else:
        return os.environ.get("HOME")


def get_swanlab_folder() -> str:
    """获取用户家目录的.swanlab文件夹路径，如果不存在此文件夹就创建

    Returns
    -------
    str
        用户家目录的.swanlab文件夹路径
    """
    user_home = get_user_home()
    swanlab_folder = os.path.join(user_home, ".swanlab")
    if not os.path.exists(swanlab_folder):
        os.mkdir(swanlab_folder)
    return swanlab_folder
