"""
@author: cunyue
@file: experiment.py
@time: 2025/12/11 18:37
@description: 定义实验相关的后端API接口
"""

from typing import Literal, Dict, TYPE_CHECKING, List, Union

from swanlab.core_python.api.type import RunType
from .utils import to_camel_case, parse_column_type

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


def send_experiment_heartbeat(
    client: "Client",
    *,
    cuid: str,
    flag_id: str,
):
    """
    发送实验心跳，保持实验处于活跃状态
    :param client: 已登录的客户端实例
    :param cuid: 实验唯一标识符
    :param flag_id: 实验标记ID
    """
    client.post(f"/house/experiments/{cuid}/heartbeat", {"flagId": flag_id})


def update_experiment_state(
    client: "Client",
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


def get_project_experiments(
    client: "Client",
    *,
    path: str,
    filters: Dict[str, object] = None,
) -> Union[List[RunType], Dict[str, List[RunType]]]:
    """
    获取指定项目下的所有实验信息
    若有实验分组，则返回一个字典，使用时需递归展平实验数据
    :param client: 已登录的客户端实例
    :param path: 项目路径 username/project
    :param filters: 筛选实验的条件，可选
    """
    parsed_filters = (
        [
            {
                "key": to_camel_case(key) if parse_column_type(key) == 'STABLE' else key.split('.', 1)[-1],
                "active": True,
                "value": [value],
                "op": 'EQ',
                "type": parse_column_type(key),
            }
            for key, value in filters.items()
        ]
        if filters
        else []
    )
    res = client.post(f"/project/{path}/runs/shows", data={'filters': parsed_filters})
    return res[0]


def get_single_experiment(client: "Client", *, path: str) -> RunType:
    """
    获取指定项目下的所有实验信息
    若有实验分组，则返回一个字典，使用时需递归展平实验数据
    :param client: 已登录的客户端实例
    :param path: 实验路径 username/project/expid
    """
    proj_path, expid = path.rsplit('/', 1)
    res = client.get(f"/project/{proj_path}/runs/{expid}")
    return res[0]


def get_experiment_metrics(client: "Client", *, expid: str, key: str) -> Dict[str, str]:
    """
    获取指定字段的指标数据，返回csv网址
    :param client: 已登录的客户端实例
    :param expid: 实验cuid
    :param key: 指定字段列表
    """
    res = client.get(f"/experiment/{expid}/column/csv", params={'key': key})
    return res[0]


__all__ = [
    "send_experiment_heartbeat",
    "update_experiment_state",
    "get_project_experiments",
    "get_single_experiment",
    "get_experiment_metrics",
]
