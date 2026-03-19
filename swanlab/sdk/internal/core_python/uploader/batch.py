"""
@author: caddiesnew
@file: batch.py
@time: 2026/3/19
@description: 分批上传函数
"""

import time
from typing import Callable, Dict, List, Literal, Optional, TypedDict, Union

from swanlab.sdk.internal.core_python.client import _get_client as get_client


class MetricDict(TypedDict):
    projectId: str
    experimentId: str
    type: str
    metrics: List[dict]
    flagId: Optional[str]


def create_data(metrics: List[dict], metrics_type: str) -> MetricDict:
    """携带上传日志的指标信息。"""
    return {
        "projectId": "",
        "experimentId": "",
        "type": metrics_type,
        "metrics": metrics,
        "flagId": None,
    }


def _generate_chunks(data: Union[MetricDict, Dict, List], per_request_len: int):
    """
    生成器：统一处理字典和列表的分片逻辑。
    yield: (分片数据, 分片长度)
    """
    if per_request_len == -1:
        if isinstance(data, list):
            yield data, len(data)
        elif isinstance(data, dict) and "metrics" in data:
            yield data, len(data["metrics"])
        else:
            yield data, 1
        return
    if isinstance(data, dict):
        metrics = data.get("metrics", [])
        curr_len = len(metrics)
        if curr_len <= per_request_len:
            yield data, curr_len
        else:
            for i in range(0, curr_len, per_request_len):
                chunk = metrics[i : i + per_request_len]
                yield {**data, "metrics": chunk}, len(chunk)

    elif isinstance(data, list):
        curr_len = len(data)
        if curr_len <= per_request_len:
            yield data, curr_len
        else:
            for i in range(0, curr_len, per_request_len):
                chunk = data[i : i + per_request_len]
                yield chunk, len(chunk)


def trace_metrics(
    url: str,
    data: Optional[Union[MetricDict, dict, list]] = None,
    method: Literal["post", "put"] = "post",
    per_request_len: int = 1000,
    upload_callback: Optional[Callable] = None,
):
    """
    分片指标上传方法。
    :param url: 上传 URL
    :param data: 要上传的数据
    :param method: HTTP 方法
    :param per_request_len: 每批上传的数量
    :param upload_callback: 上传进度回调函数
    """
    if data is None:
        return

    is_split_mode = False
    if per_request_len != -1:
        total_len = len(data.get("metrics", [])) if isinstance(data, dict) else len(data or [])
        if total_len == 0:
            return
        is_split_mode = total_len > per_request_len

    for chunk, chunk_len in _generate_chunks(data, per_request_len):
        http_client = get_client()
        resp = getattr(http_client, method)(url, chunk, retries=0)
        if resp and resp.raw.status_code == 202:
            pass  # 预留 pending 处理
        if upload_callback:
            upload_callback(chunk_len)
        if is_split_mode:
            time.sleep(1)


__all__ = [
    "MetricDict",
    "create_data",
    "trace_metrics",
]
