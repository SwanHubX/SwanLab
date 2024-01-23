#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-19 19:00:37
@File: swanlab/utils/package.py
@IDE: vscode
@Description:
    用于管理swanlab的包管理器的模块，做一些封装
"""
import json
import os


def get_package_version() -> str:
    """获取swanlab的版本号

    Returns
    -------
    str
        swanlab的版本号
    """
    try:
        # 读取package.json文件
        path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "package.json")
        with open(path, "r") as f:
            return json.load(f)["version"]
    except:
        return "unknown"
