"""
@author: caddiesnew
@file: upload.py
@time: 2026/4/22 14:01
@description: 上传相关API：conda、requirements、metadata、config、console 的上传
"""

from typing import Dict

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.upload import UploadLogMetrics, UploadMetricPayload


def upload_conda(username: str, project: str, experiment_id: str, *, content: str) -> None:
    """
    上传 conda 环境信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: conda.yaml 的原始文本内容
    """
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"conda": content})


def upload_requirements(username: str, project: str, experiment_id: str, *, content: str) -> None:
    """
    上传 Python 依赖信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: requirements.txt 的原始文本内容
    """
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"requirements": content})


def upload_metadata(username: str, project: str, experiment_id: str, *, content: Dict) -> None:
    """
    上传实验元数据。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: 元数据字典
    """
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"metadata": content})


def upload_config(username: str, project: str, experiment_id: str, *, content: Dict) -> None:
    """
    上传实验配置信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: 配置字典，格式为 {key: {value, sort, desc}}
    """
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"config": content})


def upload_console(project_id: str, experiment_id: str, *, metrics: UploadLogMetrics) -> None:
    """
    上传控制台日志信息。

    参考 Legacy upload_logs，将 ConsoleRecord 转换为日志指标格式后批量发送到 /house/metrics。

    :param project_id: 所属项目 ID
    :param experiment_id: 实验唯一标识符
    :param metrics: 控制台日志指标列表
    """
    if not metrics:
        return
    data: UploadMetricPayload = {
        "projectId": project_id,
        "experimentId": experiment_id,
        "type": "log",
        "metrics": metrics,
    }
    # retries 设置为 0 表示不重试，重试机制交给sender外层实现
    client.post("/house/metrics", data, retries=0)
