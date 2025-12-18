"""
@author: cunyue
@file: upload.py
@time: 2025/6/16 14:13
@description: 定义上传函数
"""

from functools import wraps
from typing import List

from swanlab.log import swanlog
from .batch import MetricDict, create_data, trace_metrics
from .model import ColumnModel, MediaModel, ScalarModel, FileModel, LogModel
from ..api.service import upload_to_cos
from ..client import get_client, safe_request, decode_response
from ...error import ApiError


def skip_if_empty(log_message: str):
    """
    装饰器：检查第一个参数（列表）是否为空。
    如果为空，打印 debug 日志并返回 None，不再执行原函数。
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # 获取第一个位置参数（通常是 logs, media_metrics 等列表）
            # 如果使用关键字参数调用，这里需要根据实际情况调整，但在你的代码中主要似乎是位置参数
            data_list = args[0] if args else next(iter(kwargs.values())) if kwargs else []
            # 检查是否有长度且长度为0
            if hasattr(data_list, "__len__") and len(data_list) == 0:
                return swanlog.debug(log_message)

            return func(*args, **kwargs)

        return wrapper

    return decorator


HOUSE_URL = '/house/metrics'


@safe_request
@skip_if_empty("No logs to upload.")
def upload_logs(logs: List[LogModel]):
    """
    上传日志信息
    :param logs: 日志信息集合
    """
    metrics = []
    for log in logs:
        metrics.extend([{"level": log['level'], **l} for l in log['contents']])
    if len(metrics) == 0:
        return swanlog.debug("No log metrics to upload.")
    data = create_data(metrics, "log")
    trace_metrics(HOUSE_URL, data)
    return None


@safe_request
@skip_if_empty("No media metrics to upload.")
def upload_media_metrics(media_metrics: List[MediaModel]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据集合
    """
    client = get_client()
    buffers = []
    for media in media_metrics:
        media.buffers and buffers.extend(media.buffers)
    # TODO 与 trace_metrics一样，注释掉 pending 判断
    # if not client.pending:
    if len(buffers) > 0:
        upload_to_cos(client, cuid=client.exp_id, buffers=buffers)
    # 上传指标信息
    trace_metrics(HOUSE_URL, create_data([x.to_dict() for x in media_metrics], MediaModel.type.value))


@safe_request
@skip_if_empty("No scalar metrics to upload.")
def upload_scalar_metrics(scalar_metrics: List[ScalarModel]):
    """
    上传指标的标量数据
    """
    data = create_data([x.to_dict() for x in scalar_metrics], ScalarModel.type.value)
    # 上传指标信息
    trace_metrics(HOUSE_URL, data)


@safe_request
@skip_if_empty("No files to upload.")
def upload_files(files: List[FileModel]):
    """
    上传files文件夹中的内容
    :param files: 文件列表，内部为文件信息（text/dict）
    """
    http = get_client()
    # 去重所有的FileModel，留下一个
    file_model = FileModel.create(files)
    # 如果没有文件需要上传，直接返回
    if file_model.empty:
        return
    data = file_model.to_dict()
    trace_metrics(
        f'/project/{http.groupname}/{http.projname}/runs/{http.exp_id}/profile',
        data,
        method="put",
        per_request_len=-1,
    )
    return


@safe_request
@skip_if_empty("No columns to upload.")
def upload_columns(columns: List[ColumnModel]):
    """
    批量上传并创建 columns，每个请求的列长度有一个最大值
    """
    http = get_client()
    url = f'/experiment/{http.exp_id}/columns'
    # 分批上传
    try:
        trace_metrics(url, [x.to_dict() for x in columns], per_request_len=3000)
    except ApiError as e:
        # 处理实验不存在的异常
        if e.resp.status_code == 404:
            resp = decode_response(e.resp)
            # 实验不存在，那么剩下的列也没有必要上传了，直接返回
            if isinstance(resp, dict) and resp.get('code') == 'Disabled_Resource':
                swanlog.warning(f"Experiment {http.exp_id} has been deleted, skipping column upload.")
                return
        raise e


__all__ = [
    "MetricDict",
    "upload_logs",
    "upload_media_metrics",
    "upload_scalar_metrics",
    "upload_files",
    "upload_columns",
]
