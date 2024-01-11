#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-12 01:47:34
@File: swanlab\data\utils\file.py
@IDE: vscode
@Description:
    文件/格式检查和操作
"""
import os

# ---------------------------------- 路径操作 ----------------------------------


def check_dir_and_create(path: str) -> str:
    """检查路径是否是一个绝对路径文件夹
    检查逻辑是：
    1. 检查path输入是否是一个字符串，如果不是，抛出ValueError
    2. 检查path是否是一个绝对路径，如果不是，转换为绝对路径
    3. 创建文件夹，如果文件夹已经存在，则不做任何操作
    4. 判断目标文件夹的权限，如果没有写入权限，则抛出IOError

    Parameters
    ----------
    path : str
        待检查的路径

    Returns
    -------
    str
        将path转换为绝对路径后返回

    Raises
    ------
    ValueError
        path不是一个字符串
    IOError
        path没有写权限
    """
    if not isinstance(path, str):
        raise ValueError("path must be a string")
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    # 如果创建失败，也是抛出IOError
    try:
        os.makedirs(path, exist_ok=True)
    except Exception as e:
        raise IOError(f"create path: {path} failed, error: {e}")
    if not os.access(path, os.W_OK):
        raise IOError(f"no write permission for path: {path}")
    return path
