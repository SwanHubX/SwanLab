"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API
"""

from typing import List, Literal, Optional, Tuple, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.run.v1.run_pb2 import RUN_STATE_ABORTED, RUN_STATE_CRASHED, RunState
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.typings.core_python.api.experiment import (
    InitExperimentType,
    ResumeExperimentSummaryType,
)
from swanlab.sdk.typings.run import ResumeType
from swanlab.utils.time import parse_timestamp_s


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
) -> Tuple[InitExperimentType, bool]:
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
        if run_id is None:
            raise ValueError("run_id must be provided when resume is 'must'")
        try:
            client.get(f"/project/{username}/{project}/runs/{run_id}")
        except ApiError as e:
            if e.response.status_code == 404 and e.response.reason == "Not Found":
                raise RuntimeError(f"Experiment {run_id} does not exist in project {project}")
    labels = [{"name": tag} for tag in tags] if tags else []
    created_at_for_request = created_at.ToDatetime().isoformat() + "Z"
    body = {
        "name": name,
        "description": description,
        "createdAt": created_at_for_request,
        "colors": [color, color],
        "labels": labels if len(labels) else None,
        "job": job_type,
        "cluster": group,
        "cuid": run_id,
    }
    resp = client.post(f"/project/{username}/{project}/experiment", helper.strip_none(body, strip_empty_str=True))
    # 200代表实验已存在，开启更新模式
    # 201代表实验不存在，新建实验
    # NOTE: 后端返回值没有携带createdAt字段，如果实验不存在则使用前端传入的createdAt字段作为实验创建时间，否则再请求一次获取实验创建时间
    experiment: InitExperimentType = resp.data
    is_new_experiment = resp.raw.status_code == 201
    if is_new_experiment:
        experiment["createdAt"] = created_at_for_request
    elif not experiment.get("createdAt"):
        # 旧后端 POST 响应未携带 createdAt，回退到 GET 获取
        exp_resp = client.get(f"/project/{username}/{project}/runs/{experiment['cuid']}")
        created_at_for_response = exp_resp.data.get("createdAt", "")
        if not created_at_for_response:
            raise ValueError(
                f"Backend did not return createdAt for experiment {experiment['cuid']}; please upgrade swanlab-server"
            )
        experiment["createdAt"] = created_at_for_response
    return resp.data, resp.raw.status_code == 201


def get_experiment_summary(
    project_id: str,
    experiment_id: str,
    created_at: Union[int, str],
) -> ResumeExperimentSummaryType:
    """
    获取实验摘要
    :param project_id: 所属项目ID
    :param experiment_id: 所属实验ID
    :param created_at: 实验创建时间（ISO 8601 字符串或秒/毫秒时间戳），作为数据入库时间的查询下界
    :raises ValueError: created_at 缺失或无法解析时抛出，避免回退到全历史扫描导致慢查询
    """
    ts = parse_timestamp_s(created_at)
    resp = client.get(f"/house/metrics/summaries/{project_id}/{experiment_id}", {"all": True, "createdAt": ts})
    return resp.data


def stop_experiment(username: str, project: str, experiment_id: str, *, state: RunState, finished_at: Timestamp):
    """
    停止实验
    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 所属实验名称
    :param state: 实验状态
    :param finished_at: 实验结束时间
    """
    this_state: Literal["FINISHED", "CRASHED", "ABORTED"] = "FINISHED"
    if state == RUN_STATE_CRASHED:
        this_state = "CRASHED"
    elif state == RUN_STATE_ABORTED:
        this_state = "ABORTED"
    client.put(
        f"/project/{username}/{project}/runs/{experiment_id}/state",
        {
            "state": this_state,
            "finishedAt": finished_at.ToDatetime().isoformat() + "Z",
            "from": "sdk",
        },
    )


def send_experiment_heartbeat(*, experiment_id: str):
    """
    发送实验心跳，保持实验处于活跃状态
    :param experiment_id: 实验唯一标识符
    """
    client.post(f"/house/experiments/{experiment_id}/heartbeat", {"flagId": "123456"})
