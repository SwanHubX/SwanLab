#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
"""
import os
from typing import List, MutableMapping, Optional, Union

Env = Optional[MutableMapping]

# '描述' = "key"
# ---------------------------------- 基础环境变量 ----------------------------------

ROOT = "SWANLAB_ROOT"
"""命令执行目录SWANLAB_ROOT，在这个目录下寻找swanlog文件夹"""

# ---------------------------------- 定义变量访问方法 ----------------------------------


def get_runtime_root(env: Optional[Env] = None) -> Optional[str]:
    """获取运行时根路径

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map, by default None

    Returns
    -------
    Optional[str]
        根路径
    """
    default: Optional[str] = os.getcwd()
    if env is None:
        env = os.environ
    # ROOT拿到的应该是一个绝对路径
    return env.get(ROOT, default=default)


def get_swanlog_dir() -> Optional[str]:
    """获取swanlog路径，这是一个计算变量，
    通过`get_runtime_root()`返回值得到

    Returns
    -------
    Optional[str]
        swanlog目录路径
    """
    path = os.path.join(get_runtime_root(), "swanlog")
    if not os.path.exists(path):
        os.makedirs(path)
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
    return os.path.join(get_swanlog_dir(), "project.json")
