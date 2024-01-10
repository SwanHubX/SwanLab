#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
"""
import os
from typing import MutableMapping, Optional
from .utils.file import is_port, is_ipv4, is_abs_dir

Env = Optional[MutableMapping]

_env = dict()
"""运行时环境变量参数存储，实际上就是一个字典"""

# '描述' = "key"
# ---------------------------------- 基础环境变量 ----------------------------------

ROOT = "SWANLAB_ROOT"
"""命令执行目录SWANLAB_ROOT，在这个目录下寻找swanlog文件夹"""

PORT = "SWANLAB_SERVER_PORT"
"""服务端口SWANLAB_SERVER_PORT，服务端口"""

HOST = "SWANLAB_SERVER_HOST"
"""服务端口SWANLAB_SERVER_PORT，服务地址"""


def get_runtime_root(env: Optional[Env] = None) -> Optional[str]:
    """获取运行时根路径

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map,可以是任意实现了MutableMapping的对象, 默认将使用os.environ
    Returns
    -------
    Optional[str]
        根路径
    """
    default: Optional[str] = os.getcwd()
    # 第一次调用时，从环境变量中提取，之后就不再提取，而是从缓存中提取
    if _env.get(ROOT) is not None:
        return _env.get(ROOT)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    # ROOT拿到的应该是一个绝对路径
    path = env.get(ROOT, default=default)
    _env[ROOT] = path
    # 参数检查
    if not is_abs_dir(path):
        if os.path.exists(path):
            raise ValueError('SWANLAB_ROOT must be an absolute dir, now is "{path}"'.format(path=path))
        raise ValueError('SWANLAB_ROOT must be an absolute path, now is "{path}"'.format(path=path))
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
    default: Optional[int] = 5092
    # 第一次调用时，从环境变量中提取，之后就不再提取，而是从缓存中提取
    if _env.get(PORT) is not None:
        return _env.get(PORT)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
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


# ---------------------------------- 初始化基础环境变量 ----------------------------------

# 所有的初始化函数
function_list = [
    get_runtime_root,
    get_server_port,
    get_server_host,
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

LOG_DIR = "SWANLAB_LOG_DIR"
"""日志目录SWANLAB_LOG_DIR，日志文件存放在这个目录下"""
PROJECT_PATH = "SWANLAB_PROJECT_PATH"

# ---------------------------------- 定义变量访问方法 ----------------------------------


def get_swanlog_dir() -> Optional[str]:
    """获取swanlog路径，这是一个计算变量，
    通过`get_runtime_root()`返回值得到

    Returns
    -------
    Optional[str]
        swanlog目录路径
    """
    if _env.get(LOG_DIR) is not None:
        return _env.get(LOG_DIR)
    # 否则从环境变量中提取
    path = os.path.join(get_runtime_root(), "swanlog")
    if not os.path.exists(path):
        os.makedirs(path)
    _env[LOG_DIR] = path
    return path


# TODO 后续改为数据库路径
def get_runtime_project() -> Optional[str]:
    """获取运行时项目配置，这是一个计算变量，
    通过`get_swanlog_dir()`返回值得到

    Returns
    -------
    Optional[str]
        项目配置文件路径
    """
    if _env.get(PROJECT_PATH) is not None:
        return _env.get(PROJECT_PATH)
    # 否则从环境变量中提取
    _env[PROJECT_PATH] = os.path.join(get_swanlog_dir(), "project.json")
    return _env.get(PROJECT_PATH)
