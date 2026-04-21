"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API
"""

from typing import Dict, List, Literal, Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.run.v1.run_pb2 import RUN_STATE_ABORTED, RUN_STATE_CRASHED, RunState
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.typings.core_python.api.experiment import InitExperimentType, RunType
from swanlab.sdk.typings.run import ResumeType, RunStateType
from swanlab.utils import parse_column_type, to_camel_case


def create_or_resume_experiment(
    username: str,
    project: str,
    *,
    name: str,
    resume: ResumeType,
    run_id: Optional[str] = None,
    color: str,
    description: Optional[str],
    job_type: Optional[str],
    group: Optional[str],
    tags: Optional[List[str]],
    created_at: Timestamp,
) -> InitExperimentType:
    """
    初始化实验，获取存储信息
    :param username: 所属用户名
    :param project: 所属项目名称
    :param name: 所属实验名称
    :param resume: 恢复上一次实验的状态
    :param run_id: 上一次实验的ID
    :param color: 实验颜色
    :param description: 实验描述
    :param job_type: 任务类型
    :param group: 实验组
    :param tags: 实验标签
    :param created_at: 实验创建时间，格式为 ISO 8601
    """
    if resume == "must":
        assert run_id is not None, "run_id must be provided when resume is 'must'"
        try:
            client.get(f"/project/{username}/{project}/runs/{run_id}")
        except ApiError as e:
            if e.response.status_code == 404 and e.response.reason == "Not Found":
                raise RuntimeError(f"Experiment {run_id} does not exist in project {project}")
    labels = [{"name": tag} for tag in tags] if tags else []
    body = {
        "name": name,
        "description": description,
        "createdAt": created_at.ToDatetime().isoformat() + "Z",
        "colors": [color, color],
        "labels": labels if len(labels) else None,
        "job": job_type,
        "cluster": group,
        "cuid": run_id,
    }
    resp = client.post(f"/project/{username}/{project}/experiment", helper.strip_none(body, strip_empty_str=True))
    # 200代表实验已存在，开启更新模式
    # 201代表实验不存在，新建实验
    return resp.data


def stop_experiment(username: str, project: str, cuid: str, *, state: RunState, finished_at: Timestamp):
    """
    停止实验
    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 所属实验名称
    :param state: 实验状态
    :param finished_at: 实验结束时间
    """
    this_state: Literal["FINISHED", "CRASHED", "ABORTED"] = "FINISHED"
    if state == RUN_STATE_CRASHED:
        this_state = "CRASHED"
    elif state == RUN_STATE_ABORTED:
        this_state = "ABORTED"
    resp = client.put(
        f"/project/{username}/{project}/runs/{cuid}/state",
        {
            "state": this_state,
            "finishedAt": finished_at.ToDatetime().isoformat() + "Z",
            "from": "sdk",
        },
    )
    return resp.raw.status_code == 201


def send_experiment_heartbeat(*, cuid: str, flag_id: str) -> None:
    """
    发送实验心跳，保持实验处于活跃状态
    :param cuid: 实验唯一标识符
    :param flag_id: 实验标记ID
    """
    client.post(f"/house/experiments/{cuid}/heartbeat", {"flagId": flag_id})


def update_experiment_state(
    *,
    username: str,
    projname: str,
    cuid: str,
    state: RunStateType,
    finished_at: Optional[str] = None,
) -> None:
    """
    更新实验状态
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
    put_data = {k: v for k, v in put_data.items() if v is not None}
    client.put(f"/project/{username}/{projname}/runs/{cuid}/state", put_data)


def get_project_experiments(
    *,
    path: str,
    filters: Optional[Dict[str, object]] = None,
) -> Union[List[RunType], Dict[str, List[RunType]]]:
    """
    获取指定项目下的所有实验信息
    若有实验分组，则返回一个字典，使用时需递归展平实验数据
    :param path: 项目路径 username/project
    :param filters: 筛选实验的条件，可选
    """
    parsed_filters = (
        [
            {
                "key": to_camel_case(key) if parse_column_type(key) == "STABLE" else key.split(".", 1)[-1],
                "active": True,
                "value": [value],
                "op": "EQ",
                "type": parse_column_type(key),
            }
            for key, value in filters.items()
        ]
        if filters
        else []
    )
    return client.post(f"/project/{path}/runs/shows", data={"filters": parsed_filters}).data


def get_single_experiment(*, path: str) -> RunType:
    """
    获取指定实验信息
    :param path: 实验路径 username/project/expid
    """
    proj_path, expid = path.rsplit("/", 1)
    return client.get(f"/project/{proj_path}/runs/{expid}").data


def get_experiment_metrics(*, expid: str, key: str) -> Dict[str, str]:
    """
    获取指定字段的指标数据，返回csv网址
    :param expid: 实验cuid
    :param key: 指定字段列表
    """
    return client.get(f"/experiment/{expid}/column/csv", params={"key": key}).data


def delete_experiment(*, path: str) -> None:
    """
    删除指定实验
    :param path: 实验路径 'username/project/expid'
    """
    proj_path, expid = path.rsplit("/", 1)
    client.delete(f"/project/{proj_path}/runs/{expid}")
