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
import enum
import sys

Env = Optional[MutableMapping]

_env = dict()
"""运行时环境变量参数存储，实际上就是一个字典"""

# ---------------------------------- 基础环境变量 ----------------------------------
# '描述' = "key"
ROOT = "SWANLAB_LOG_DIR"
"""命令执行目录SWANLAB_LOG_DIR，日志文件存放在这个目录下，如果自动生成，则最后的目录名为swanlog，第一次调用时如果路径不存在，会自动创建路径"""

PORT = "SWANLAB_SERVER_PORT"
"""服务端口SWANLAB_SERVER_PORT，服务端口"""

HOST = "SWANLAB_SERVER_HOST"
"""服务端口SWANLAB_SERVER_PORT，服务地址"""

DEV = "SWANLAB_DEV"
"""是否是开发模式SWANLAB_DEV，开发模式下会打印更多的日志信息，并且切换一些配置"""

MODE = "SWANLAB_MODE"
"""运行模式SWANLAB_MODE，有cloud、local、disabled等模式，针对Callback回调的不同实例化方式，具体效果为：
1. disabled: 禁用模式，不会上传日志到云端，也不会记录日志，swanlab只会执行解析log的功能，不会输出任何内容到本地
2. cloud: 云端模式，是云端与本地的混合模式，会上传日志到云端，也会记录日志到本地
3. cloud-only: 云端模式，只会上传日志到云端，不会记录日志到本地，但是会生成一些临时文件
4. local: 本地模式，不会上传日志到云端，使用swanlab本地版本
此外，SWANLAB_MODE为disabled时，会开启非严格模式，非严格模式不再要求文件路径存在
"""


class SwanLabMode(enum.Enum):
    DISABLED = "disabled"
    CLOUD = "cloud"
    # CLOUD_ONLY = "cloud-only"
    LOCAL = "local"


def get_mode(env: Optional[Env] = None) -> Optional[str]:
    """
    获取运行模式，返回值为SwanLabMode枚举value
    """
    if _env.get(MODE) is not None:
        return _env.get(MODE)
    # 否则从环境变量中提取
    if env is None:
        env = os.environ
    default: Optional[str] = "cloud"
    mode = env.get(MODE, default)
    allowed = [mode.value for mode in SwanLabMode]
    if mode not in allowed:
        raise ValueError('SWANLAB_MODE must be one of {allowed}, now is "{mode}"'.format(allowed=allowed, mode=mode))
    _env[MODE] = mode
    return mode


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
    # 严格模式路径必须存在
    assert_exist(
        path,
        target_type="folder",
        desc=(
            'The log file was not found in the default path "{path}". '
            'Please use the "swanlab watch -l <LOG '
            'PATH>" command to specify the location of the log path."'.format(path=path)
            if path == default
            else 'SWANLAB_LOG_DIR must be an existing path, now is "{path}"'.format(path=path)
        ),
        t_desc='SWANLAB_LOG_DIR must be a directory, now is "{path}"'.format(path=path),
    )
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
function_list = [get_mode, get_swanlog_dir, get_server_port, get_server_host, is_dev]


def init_env(env: Optional[Env] = None):
    """初始化环境变量

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map,可以是任意实现了MutableMapping的对象, 默认将使用os.environ
    """
    reset_env()
    for func in function_list:
        func(env)


def reset_env():
    """重置"""
    _env.clear()


# ---------------------------------- 定义计算变量访问方法 ----------------------------------


def is_strict_mode() -> bool:
    """
    是否是严格模式，严格模式下会要求文件路径存在，否则会抛出异常
    """
    return get_mode() != SwanLabMode.DISABLED.value


def get_db_path() -> Optional[str]:
    """
    获取数据库路径，这是一个计算变量，每次调用都会重新计算
    """
    return os.path.join(get_swanlog_dir(), "runs.swanlab")


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
    try:
        if not assert_exist(swanlab_folder, ra=False, target_type="folder"):
            os.mkdir(swanlab_folder)
    except NotADirectoryError:
        os.remove(swanlab_folder)
        os.mkdir(swanlab_folder)
    return swanlab_folder


def assert_exist(path: str, target_type: str = None, ra: bool = True, desc: str = None, t_desc: str = None) -> bool:
    """
    检查文件是否存在，严格模式下，文件不存在会抛出异常，或者可以手动通过参数控制，存在则返回True，否则返回False
    :param path: 文件路径
    :param target_type: 文件类型(folder, file)，如果文件类型与预期不符，会抛出异常，非严格模式下不检测，为None不检测文件类型
    :param ra: 文件不存在时是否抛出异常，非严格模式下强制不抛出，此参数无效
    :param desc: 异常描述信息，非严格模式下强制不抛出，此参数无效
    :param t_desc: 文件类型描述信息，非严格模式下强制不抛出，此参数无效
    """
    if not is_strict_mode():
        return os.path.exists(path)
    if not os.path.exists(path):
        if ra:
            raise FileNotFoundError(desc or "{path} not existed".format(path=path))
        else:
            return False
    # 检查文件类型
    if target_type is not None:
        if target_type == "folder":
            if not os.path.isdir(path):
                raise NotADirectoryError(t_desc or "{path} is not a folder".format(path=path))
        elif target_type == "file":
            if not os.path.isfile(path):
                raise IsADirectoryError(t_desc or "{path} is not a file".format(path=path))
    return True
