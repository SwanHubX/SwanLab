"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API
"""

from datetime import datetime, timezone
from typing import List, Literal, Optional, Tuple, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.run.v1.run_pb2 import RUN_STATE_ABORTED, RUN_STATE_CRASHED, RunState
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.typings.core_python.api.experiment import (
    CREATED_AT_MAX,
    CREATED_AT_MIN,
    InitExperimentType,
    ResumeExperimentSummaryType,
)
from swanlab.sdk.typings.run import ResumeType


def parse_timestamp_s(value: Union[int, str, None]) -> int:
    """将时间统一转换为秒级 Unix 时间戳（10 位），用于 House 查询接口的 ``createdAt`` 参数。

    支持格式：
    - int / 数字字符串：秒（10 位）或毫秒（13 位）级时间戳，自动归一化为秒
    - ISO 8601 字符串：如 ``"2024-08-01T00:00:00Z"``、``"2024-08-01T08:00:00+08:00"``，无时区时按 UTC 解析
    - None / 空字符串 / 无法解析：返回 0
    """
    if value is None:
        return 0
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return 0
        if not text.isdigit():
            try:
                dt = datetime.fromisoformat(text.replace("Z", "+00:00"))
            except ValueError:
                return 0
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            return int(dt.timestamp())
        value = int(text)
    if not isinstance(value, int) or value <= 0:
        return 0
    while value > 9_999_999_999:
        value //= 1000
    return value


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
    return resp.data, resp.raw.status_code == 201


def get_experiment_summary(
    project_id: str,
    experiment_id: str,
    created_at: Union[int, str, None] = None,
) -> ResumeExperimentSummaryType:
    """
    获取实验摘要
    :param project_id: 所属项目ID
    :param experiment_id: 所属实验ID
    :param created_at: 实验创建时间（ISO 8601 字符串或秒/毫秒时间戳），作为数据入库时间的查询下界；
        无法解析时回退到下界值（全历史扫描）以保证查询正确性
    """
    ts = parse_timestamp_s(created_at)
    if not CREATED_AT_MIN <= ts <= CREATED_AT_MAX:
        ts = CREATED_AT_MIN
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
