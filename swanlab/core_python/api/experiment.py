"""
@author: cunyue
@file: experiment.py
@time: 2025/12/11 18:37
@description: 定义实验相关的后端API接口
"""

from typing import Literal

from swanlab.core_python.client import Client


def update_experiment_state(
    client: Client,
    *,
    username: str,
    projname: str,
    cuid: str,
    state: Literal['FINISHED', 'CRASHED', 'ABORTED'],
    finished_at: str = None,
):
    """
    更新实验状态，注意此接口会将客户端标记为 pending 状态，表示实验已结束
    :param client: 已登录的客户端实例
    :param username: 实验所属用户名
    :param projname: 实验所属项目名称
    :param cuid: 实验唯一标识符
    :param state: 实验状态
    :param finished_at: 实验结束时间，格式为 ISO 8601，如果不提供则使用当前时间
    """
    put_data = {
        "state": state,
        "finishedAt": finished_at,
        "from": "sdk",
    }
    put_data = {k: v for k, v in put_data.items() if v is not None}  # 移除值为None的键
    client.put(f"/project/{username}/{projname}/runs/{cuid}/state", put_data)
    client.pending = True
