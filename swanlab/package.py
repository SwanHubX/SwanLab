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
from .env import is_dev, get_swanlab_folder, is_strict_mode
from .utils.key import get_key
from .error import KeyFileError

package_path = None
if is_dev():
    # 当前运行时的位置
    package_path = os.environ['SWANLAB_PACKAGE_PATH']
else:
    package_path = os.path.join(os.path.dirname(__file__), "package.json")


def get_package_version(p=package_path) -> str:
    """获取swanlab的版本号

    Parameters
    ----------
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        swanlab的版本号
    """
    # 读取package.json文件
    with open(p, "r") as f:
        return json.load(f)["version"]


def get_host_web(p: str = package_path) -> str:
    """获取swanlab网站网址

    Parameters
    ----------
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        swanlab网站的网址
    """
    with open(p, "r") as f:
        return json.load(f)["host"]["web"]


def get_host_api(p: str = package_path) -> str:
    """获取swanlab网站api网址

    Parameters
    ----------
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        swanlab网站的api网址
    """
    with open(p, "r") as f:
        return json.load(f)["host"]["api"]


def get_user_setting_path(p: str = package_path) -> str:
    """获取用户设置的url

    Parameters
    ----------
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        用户设置的url
    """
    return get_host_web(p) + "/settings"


def get_project_url(username: str, projname: str, p: str = package_path) -> str:
    """获取项目的url

    Parameters
    ----------
    username : str
        用户名
    projname : str
        项目名
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        项目的url
    """
    return get_host_web(p) + "/" + username + "/" + projname


def get_experiment_url(username: str, projname: str, expid: str, p: str = package_path) -> str:
    """获取实验的url

    Parameters
    ----------
    username : str
        用户名
    projname : str
        项目名
    expid : str
        实验id
    p : str, optional
        package.json文件路径，默认为项目的package.json

    Returns
    -------
    str
        实验的url
    """
    return get_project_url(username, projname, p) + "/" + expid


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
    if not is_strict_mode():
        return
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
                    info = ("The version of logdir is old (Created by swanlab<=0.1.4), the current version of "
                            "SwanLab doesn't support this logfile. If you need to watch this logfile, please use the "
                            "transfer script: https://github.com/SwanHubX/SwanLab/blob/main/script/transfer_logfile_0"
                            ".1.4.py'")
                elif mode == "init":
                    info = ("The version of logdir is old (Created by swanlab<=0.1.4), the current version of "
                            "SwanLab doesn't support this logfile. If you need to continue train in this logdir, "
                            "please use the transfer script: "
                            "https://github.com/SwanHubX/SwanLab/blob/main/script/transfer_logfile_0.1.4.py'")
                else:
                    info = "version_limit function only support mode in ['watch', 'init']"
                raise ValueError(info)
    return


def is_login() -> bool:
    """判断是否已经登录

    Returns
    -------
    bool
        是否已经登录，但是不保证登录的key可用
    """
    try:
        _ = get_key(os.path.join(get_swanlab_folder(), ".netrc"), get_host_api())[2]
        return True
    except KeyFileError:
        return False
