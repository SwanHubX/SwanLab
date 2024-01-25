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


def version_limit(path: str, mode: str) -> None:
    """限制版本号在v0.1.5之后
    主要原因是因为在v0.1.5之后使用了数据库连接，而在之前的版本中并没有使用数据库连接
    本函数用于检查swanlog文件夹内容是否是v0.1.5之后的版本

    Parameters
    ----------
    path : str
        swanlog目录
    mode : str
        模式，只能是['watch', 'init']中的一个，分别对应swanlab watch和swanlab init
        主要是显示不同的错误信息

    Raises
    ------
    ValueError
        如果版本号低于v0.1.5则抛出异常
    """
    # 判断文件夹内是否存在runs.swanlog文件
    if os.path.exists(os.path.join(path, "runs.swanlog")):
        return
    # 不存在则进一步判断，如果存在project.json文件，且可以被json读取，并且存在version字段，则说明版本号小于v0.1.5
    if os.path.exists(os.path.join(path, "project.json")):
        with open(os.path.join(path, "project.json"), "r") as f:
            project = json.load(f)
            # print(project.get("version"))
            if project.get("version") is not None:
                # 报错，当前目录只允许v0.1.5之前的版本，请降级到v0.1.4
                if mode == "watch":
                    info = "The version of logdir's file is old (Created by swanlab<=0.1.4), the current version of SwanLab doesn't support this logfile. If you need to watch this logfile, please use the transfer script: https://github.com/SwanHubX/SwanLab/blob/main/script/transfer_logfile_0.1.4.py'"
                elif mode == "init":
                    info = "The version of logdir's file is old (Created by swanlab<=0.1.4), the current version of SwanLab doesn't support this logfile. If you need to continue train in this logfir, please use the transfer script: https://github.com/SwanHubX/SwanLab/blob/main/script/transfer_logfile_0.1.4.py'"
                else:
                    info = "version_limit function only support mode in ['watch', 'init']"
                raise ValueError(info)
    return
