"""
@author: cunyue
@file: upload.py
@time: 2025/6/16 14:13
@description: 定义上传函数
"""

from typing import List

from swanlab.log import swanlog
from .batch import MetricDict, create_data, trace_metrics
from .model import ColumnModel, MediaModel, ScalarModel, FileModel, LogModel
from ..api.service import upload_to_cos
from ..client import get_client, safe_request, decode_response
from ...error import ApiError

HOUSE_URL = '/house/metrics'


@safe_request
def upload_logs(logs: List[LogModel]):
    """
    上传日志信息
    :param logs: 日志信息集合
    """
    metrics = []
    for log in logs:
        metrics.extend([{"level": log['level'], **l} for l in log['contents']])
    if len(metrics) == 0:
        return swanlog.debug("No logs to upload.")
    data = create_data(metrics, "log")
    trace_metrics(HOUSE_URL, data)
    return None


@safe_request
def upload_media_metrics(media_metrics: List[MediaModel]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据集合
    """
    client = get_client()
    buffers = []
    for media in media_metrics:
        media.buffers and buffers.extend(media.buffers)
    if not client.pending:
        upload_to_cos(client, cuid=client.exp_id, buffers=buffers)
        # 上传指标信息
        trace_metrics(HOUSE_URL, create_data([x.to_dict() for x in media_metrics], MediaModel.type.value))


@safe_request
def upload_scalar_metrics(scalar_metrics: List[ScalarModel]):
    """
    上传指标的标量数据
    """
    data = create_data([x.to_dict() for x in scalar_metrics], ScalarModel.type.value)
    trace_metrics(HOUSE_URL, data)


@safe_request
def upload_files(files: List[FileModel]):
    """
    上传files文件夹中的内容
    :param files: 文件列表，内部为文件信息（text/dict）
    """
    http = get_client()
    # 去重所有的FileModel，留下一个
    if len(files) == 0:
        swanlog.warning("No files to upload.")
        return
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
def upload_columns(columns: List[ColumnModel]):
    """
    批量上传并创建 columns，每个请求的列长度有一个最大值
    """
    http = get_client()
    url = f'/experiment/{http.exp_id}/columns'
    # 如果列表长度为0，则跳过
    if len(columns) == 0:
        swanlog.debug("No columns to upload.")
        return
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
