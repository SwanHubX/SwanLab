"""
@author: caddiesnew
@file: environment.py
@time: 2026/4/22 14:01
@description: 环境相关API：conda、requirements、metadata、config 的上传
"""

from typing import Dict

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper


def upload_conda(username: str, project: str, cuid: str, *, content: str) -> None:
    """
    上传 conda 环境信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 实验唯一标识符
    :param content: conda.yaml 的原始文本内容
    """
    client.put(
        f"/project/{username}/{project}/runs/{cuid}/profile",
        helper.strip_none({"conda": content}),
    )


def upload_requirements(username: str, project: str, cuid: str, *, content: str) -> None:
    """
    上传 Python 依赖信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 实验唯一标识符
    :param content: requirements.txt 的原始文本内容
    """
    client.put(
        f"/project/{username}/{project}/runs/{cuid}/profile",
        helper.strip_none({"requirements": content}),
    )


def upload_metadata(username: str, project: str, cuid: str, *, content: Dict) -> None:
    """
    上传实验元数据。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 实验唯一标识符
    :param content: 元数据字典
    """
    client.put(
        f"/project/{username}/{project}/runs/{cuid}/profile",
        helper.strip_none({"metadata": content}),
    )


def upload_config(username: str, project: str, cuid: str, *, content: Dict) -> None:
    """
    上传实验配置信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 实验唯一标识符
    :param content: 配置字典，格式为 {key: {value, sort, desc}}
    """
    client.put(
        f"/project/{username}/{project}/runs/{cuid}/profile",
        helper.strip_none({"config": content}),
    )
