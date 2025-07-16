"""
@author: cunyue
@file: upload.py
@time: 2025/6/16 14:13
@description: 定义上传函数
"""

from typing import List, Union, Literal

from swanlab.log import swanlog
from .model import ColumnModel, MediaModel, ScalarModel, FileModel, LogModel
from ..client import get_client, sync_error_handler, decode_response
from ...error import ApiError

house_url = '/house/metrics'


def create_data(metrics: List[dict], metrics_type: str) -> dict:
    """
    携带上传日志的指标信息
    """
    client = get_client()
    # Move 等实验需要将数据上传到根实验上
    exp_id = client.exp.root_exp_cuid or client.exp.cuid
    proj_id = client.exp.root_proj_cuid or client.proj.cuid

    flag_id = client.exp.flag_id
    return {
        "projectId": proj_id,
        "experimentId": exp_id,
        "type": metrics_type,
        "metrics": metrics,
        "flagId": flag_id,
    }


def trace_metrics(url: str, data: Union[dict, list] = None, method: Literal['post', 'put'] = 'post'):
    """
    创建指标数据方法，如果 client 处于挂起状态，则不进行上传
    :param url: 上传的URL地址
    :param data: 上传的数据，可以是字典或列表
    :param method: 请求方法，默认为 'post'
    """
    client = get_client()
    if client.pending:
        return
    _, resp = getattr(client, method)(url, data)
    if resp.status_code == 202:
        client.pending = True


@sync_error_handler
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
    trace_metrics(house_url, data)
    return None


@sync_error_handler
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
        client.upload_files(buffers)
        # 上传指标信息
        trace_metrics(house_url, create_data([x.to_dict() for x in media_metrics], MediaModel.type.value))


@sync_error_handler
def upload_scalar_metrics(scalar_metrics: List[ScalarModel]):
    """
    上传指标的标量数据
    """
    data = create_data([x.to_dict() for x in scalar_metrics], ScalarModel.type.value)
    trace_metrics(house_url, data)


@sync_error_handler
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
    trace_metrics(f'/project/{http.groupname}/{http.projname}/runs/{http.exp_id}/profile', data, method="put")
    return


@sync_error_handler
def upload_columns(columns: List[ColumnModel], per_request_len: int = 3000):
    """
    批量上传并创建 columns，每个请求的列长度有一个最大值
    """
    http = get_client()
    url = f'/experiment/{http.exp_id}/columns'
    # 将columns拆分成多个小的列表，每个列表的长度不能超过单个请求的最大长度
    columns_list = []
    columns_count = len(columns)
    for i in range(0, columns_count, per_request_len):
        columns_list.append([columns[i + j].to_dict() for j in range(min(per_request_len, columns_count - i))])
    # 上传每个列表
    for columns in columns_list:
        # 如果列表长度为0，则跳过
        if len(columns) == 0:
            continue
        try:
            trace_metrics(url, columns)
        except ApiError as e:
            # 处理实验不存在的异常
            if e.resp.status_code == 404:
                resp = decode_response(e.resp)
                # 实验不存在，那么剩下的列也没有必要上传了，直接返回
                if isinstance(resp, dict) and resp.get('code') == 'Disabled_Resource':
                    swanlog.warning(f"Experiment {http.exp_id} has been deleted, skipping column upload.")
                    return
            raise e
    return


__all__ = [
    "upload_logs",
    "upload_media_metrics",
    "upload_scalar_metrics",
    "upload_files",
    "upload_columns",
]
