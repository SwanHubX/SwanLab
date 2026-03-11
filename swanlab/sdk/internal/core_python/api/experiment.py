"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API
"""

from typing import List, Optional

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.run import ResumeType
from swanlab.sdk.utils.helper import strip_none


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
    created_at: Optional[str] = None,
) -> bool:
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
    if resume == "never":
        assert created_at is None, "created_at must be None when resume is 'never'"
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
        "createdAt": created_at,
        "colors": [color, color],
        "labels": labels if len(labels) else None,
        "job": job_type,
        "cluster": group,
        "cuid": run_id,
    }
    resp = client.post(f"/project/{username}/{project}/experiment", strip_none(body))
    # 200代表实验已存在，开启更新模式
    # 201代表实验不存在，新建实验
    return resp.raw.status_code == 201
