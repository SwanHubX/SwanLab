"""
@author: caddiesnew
@file: upload.py
@time: 2026/4/22 14:01
@description: 上传相关API：conda、requirements、metadata、config、console 的上传
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import IO, TYPE_CHECKING, Dict, List, Optional, Union

from requests.sessions import Session

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper, get_buffer_size
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.typings.core_python.api.upload import (
    UploadColumns,
    UploadLogBatch,
    UploadMediaBatch,
    UploadMetricPayload,
    UploadResource,
    UploadScalarBatch,
)

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


def upload_conda(username: str, project: str, experiment_id: str, *, content: str) -> None:
    """
    上传 conda 环境信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: conda.yaml 的原始文本内容
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"conda": content}, retries=0)


def upload_requirements(username: str, project: str, experiment_id: str, *, content: str) -> None:
    """
    上传 Python 依赖信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: requirements.txt 的原始文本内容
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"requirements": content}, retries=0)


def upload_metadata(username: str, project: str, experiment_id: str, *, content: Dict) -> None:
    """
    上传实验元数据。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: 元数据字典
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"metadata": content}, retries=0)


def upload_config(username: str, project: str, experiment_id: str, *, content: Dict) -> None:
    """
    上传实验配置信息。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param experiment_id: 实验唯一标识符
    :param content: 配置字典，格式为 {key: {value, sort, desc}}
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.put(f"/project/{username}/{project}/runs/{experiment_id}/profile", {"config": content}, retries=0)


def upload_columns(username: str, project: str, *, columns: UploadColumns) -> None:
    """
    上传列信息
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.post(f"/projects/{username}/{project}/series", columns, retries=0)


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
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
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
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
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
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    client.post("/house/metrics", data, retries=0)


def upload_resource(
    session: Session,
    experiment_id: str,
    *,
    paths: List[str],
    buffers: List[Union[IO[bytes], str, Path]],
    content_types: Optional[List[str]] = None,
    tracker: Optional[UploadTracker] = None,
) -> Dict[int, BaseException]:
    """
    上传资源文件到对象存储，先向后端请求上传凭据，然后使用凭据上传文件到对象存储。

    :param content_types: 可选的逐文件 Content-Type 列表，不传时使用二进制默认值。
    :param tracker: 可选的进度追踪器，非 None 时汇报字节级进度。
    :return: {索引: 异常} 字典，索引与 paths 对齐；全部成功时为空字典。
    """
    # retries 设置为 0 表示不重试，重试机制交给 transport 层的退避调度
    resp = client.post(
        "/resources/presigned/put",
        {"experimentId": experiment_id, "paths": paths},
        retries=0,
    )
    data = resp.data
    urls = data.get("urls") if isinstance(data, dict) else None
    if not isinstance(urls, list) or not all(isinstance(url, str) and url for url in urls):
        # 字段缺失或类型错误属于协议错误，重试也不会恢复，交由 transport 终结记录
        raise ApiError(
            code="invalid_presigned_response",
            message="presigned response missing a valid 'urls' list",
            retryable=False,
        )
    if len(urls) != len(paths):
        # 数量不匹配属于协议错误，重试也不会恢复，交由 transport 终结记录
        raise ApiError(
            code="invalid_presigned_response",
            message=f"presigned response returned {len(urls)} urls for {len(paths)} paths",
            retryable=False,
        )

    def upload_one(url: str, buffer: Union[IO[bytes], str, Path], file_key: str, size: int, content_type: str):
        # 上传单个文件，支持文件路径和内存 buffer 两种形式
        if isinstance(buffer, (str, Path)):
            with open(buffer, "rb") as f:
                _put_with_progress(session, url, f, file_key, size, tracker, content_type=content_type)
        else:
            _put_with_progress(session, url, buffer, file_key, size, tracker, content_type=content_type)

    failures: Dict[int, BaseException] = {}
    with SafeThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for index, url in enumerate(urls):
            buffer = buffers[index]
            size = get_buffer_size(buffer)
            file_key = f"{paths[index]}:{size}"
            content_type = "application/octet-stream"
            if content_types is not None and index < len(content_types):
                content_type = content_types[index]
            futures.append((executor.submit(upload_one, url, buffer, file_key, size, content_type), index, file_key))
        for future, index, file_key in futures:
            try:
                future.result()
            except Exception as e:
                _reset_tracked_file(tracker, file_key)
                failures[index] = e
    return failures


def upload_saves(
    session: Session,
    *,
    resources: Sequence[UploadResource],
    tracker: Optional[UploadTracker] = None,
) -> Dict[int, BaseException]:
    """通过预签名 URL 批量上传本地文件到对象存储，返回 {索引: 异常} 字典。"""

    def upload_one(resource: UploadResource) -> None:
        file_key = resource.get("tracker_key")
        size = resource.get("size")
        with open(resource["source_path"], "rb") as f:
            _put_with_progress(
                session, resource["url"], f, file_key, size, tracker, content_type=resource["content_type"]
            )

    failures: Dict[int, BaseException] = {}
    with SafeThreadPoolExecutor(max_workers=8) as executor:
        futures = [(executor.submit(upload_one, resource), index) for index, resource in enumerate(resources)]
        for future, index in futures:
            try:
                future.result()
            except Exception as e:
                _reset_tracked_file(tracker, resources[index].get("tracker_key"))
                failures[index] = e
    return failures


def _reset_tracked_file(tracker: Optional[UploadTracker], file_key: Optional[str]) -> None:
    """上传失败后重置文件字节级进度。"""
    if tracker is not None and file_key is not None:
        tracker.update_file_progress(file_key, file_key, 0)


def _put_with_progress(
    session: Session,
    url: str,
    file_obj: IO[bytes],
    file_key: Optional[str],
    size: Optional[int],
    tracker: Optional[UploadTracker],
    *,
    content_type: str = "application/octet-stream",
) -> None:
    """PUT 上传单个文件，有 tracker 时用 ProgressFileWrapper 汇报读字节进度。"""
    if tracker is not None and file_key is not None and size:
        data: Union[IO[bytes], ProgressFileWrapper] = ProgressFileWrapper(
            file_obj,
            on_read=lambda current: tracker.update_file_progress(file_key, file_key, current),
            total_size=size,
        )
    else:
        data = file_obj
    resp = session.put(url, data=data, headers={"Content-Type": content_type})
    resp.raise_for_status()
