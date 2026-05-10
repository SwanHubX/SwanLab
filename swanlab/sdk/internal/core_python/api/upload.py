"""
@author: caddiesnew
@file: upload.py
@time: 2026/4/22 14:01
@description: 上传相关API：conda、requirements、metadata、config、console 的上传
"""

from _io import BytesIO
from collections.abc import Sequence
from typing import Dict, List

from requests.sessions import Session

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.core_python.api.upload import (
    UploadColumns,
    UploadLogBatch,
    UploadMediaBatch,
    UploadMetricPayload,
    UploadResource,
    UploadScalarBatch,
)


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


def upload_columns(experiment_id: str, *, columns: UploadColumns) -> None:
    """
    上传列信息
    """
    client.post(f"/experiment/{experiment_id}/columns", columns, retries=0)


def upload_log(project_id: str, experiment_id: str, *, metrics: UploadLogBatch) -> None:
    """
    批量发送日志信息到 /house/metrics。

    :param project_id: 所属项目 ID
    :param experiment_id: 实验唯一标识符
    :param metrics: 控制台日志指标列表
    """
    if not metrics:
        return
    data: UploadMetricPayload = {
        "type": "log",
        "projectId": project_id,
        "experimentId": experiment_id,
        "metrics": metrics,
    }
    # retries 设置为 0 表示不重试，重试机制交给sender外层实现
    client.post("/house/metrics", data, retries=0)


def upload_scalar(project_id: str, experiment_id: str, *, metrics: UploadScalarBatch) -> None:
    """
    批量发送指标信息到 /house/metrics。

    :param project_id: 所属项目 ID
    :param experiment_id: 实验唯一标识符
    :param metrics: 控制台日志指标列表
    """
    if not metrics:
        return
    data: UploadMetricPayload = {
        "type": "scalar",
        "projectId": project_id,
        "experimentId": experiment_id,
        "metrics": metrics,
    }
    # retries 设置为 0 表示不重试，重试机制交给sender外层实现
    client.post("/house/metrics", data, retries=0)


def upload_media(project_id: str, experiment_id: str, *, metrics: UploadMediaBatch) -> None:
    """
    批量发送媒体信息到 /house/metrics。

    :param project_id: 所属项目 ID
    :param experiment_id: 实验唯一标识符
    :param metrics: 控制台日志指标列表
    """
    if not metrics:
        return
    data: UploadMetricPayload = {
        "type": "media",
        "projectId": project_id,
        "experimentId": experiment_id,
        "metrics": metrics,
    }
    # retries 设置为 0 表示不重试，重试机制交给sender外层实现
    client.post("/house/metrics", data, retries=0)


def upload_resource(session: Session, experiment_id: str, *, paths: List[str], buffers: List[BytesIO]):
    """
    上传资源文件到对象存储，先向后端请求上传凭据，然后使用凭据上传文件到对象存储。
    """
    resp = client.post(
        "/resources/presigned/put",
        {"experimentId": experiment_id, "paths": paths},
    )
    urls = resp.data["urls"]
    with SafeThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for index, url in enumerate(urls):
            futures.append(
                executor.submit(
                    session.put,
                    url,
                    data=buffers[index],
                    headers={"Content-Type": "application/octet-stream"},
                )
            )
        for future in futures:
            resp = future.result()
            resp.raise_for_status()


def upload_saves(
    session: Session,
    *,
    resources: Sequence[UploadResource],
) -> set[str]:
    """通过预签名 URL 批量上传本地文件到对象存储，返回成功上传的 source_path 集合。"""

    @safe.decorator(message="Failed to upload save file, skipping")
    def upload_one(resource: UploadResource) -> bool:
        with open(resource["source_path"], "rb") as f:
            resp = session.put(resource["url"], data=f, headers={"Content-Type": resource["content_type"]})
            resp.raise_for_status()
        return True

    uploaded: set[str] = set()
    with SafeThreadPoolExecutor(max_workers=8) as executor:
        futures = [(executor.submit(upload_one, resource), resource["source_path"]) for resource in resources]
        for future, source_path in futures:
            if future.result() is True:
                uploaded.add(source_path)

    return uploaded
