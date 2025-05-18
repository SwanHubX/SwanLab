#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:56
@File: __init__.py
@IDE: pycharm
@Description:
    上传相关接口
"""
from typing import List

from swanlab.log import swanlog
from .model import ColumnModel, MediaModel, ScalarModel, FileModel, LogModel
from ..http import get_http, sync_error_handler, decode_response
from ...error import ApiError

house_url = '/house/metrics'


def create_data(metrics: List[dict], metrics_type: str) -> dict:
    """
    携带上传日志的指标信息
    """
    http = get_http()
    return {"projectId": http.proj_id, "experimentId": http.exp_id, "type": metrics_type, "metrics": metrics}


@sync_error_handler
def upload_logs(logs: List[LogModel]):
    """
    上传日志信息
    :param logs: 日志信息集合
    """
    http = get_http()
    metrics = []
    for log in logs:
        metrics.extend([{"level": log['level'], **l} for l in log['contents']])
    if len(metrics) == 0:
        return swanlog.debug("No logs to upload.")
    data = create_data(metrics, "log")
    http.post(house_url, data)


@sync_error_handler
def upload_media_metrics(media_metrics: List[MediaModel]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据集合
    """
    http = get_http()
    buffers = []
    for media in media_metrics:
        media.buffers and buffers.extend(media.buffers)
    http.upload_files(buffers)
    # 上传指标信息
    http.post(house_url, create_data([x.to_dict() for x in media_metrics], MediaModel.type.value))


@sync_error_handler
def upload_scalar_metrics(scalar_metrics: List[ScalarModel]):
    """
    上传指标的标量数据
    """
    http = get_http()
    data = create_data([x.to_dict() for x in scalar_metrics], ScalarModel.type.value)
    http.post(house_url, data)


@sync_error_handler
def upload_files(files: List[FileModel]):
    """
    上传files文件夹中的内容
    :param files: 文件列表，内部为文件绝对路径
    """
    http = get_http()
    # 去重所有的FileModel，留下一个
    if len(files) == 0:
        return swanlog.warning("No files to upload.")
    file_model = FileModel.create(files)
    # 如果没有文件需要上传，直接返回
    if file_model.empty:
        return None
    data = file_model.to_dict()
    http.put(f'/project/{http.groupname}/{http.projname}/runs/{http.exp_id}/profile', data)
    return None


@sync_error_handler
def upload_columns(columns: List[ColumnModel], per_request_len: int = 3000):
    """
    批量上传并创建 columns，每个请求的列长度有一个最大值
    """
    http = get_http()
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
            http.post(url, columns)
        except ApiError as e:
            # 处理实验不存在的异常
            if e.resp.status_code == 404:
                resp = decode_response(e.resp)
                # 实验不存在，那么剩下的列也没有必要上传了，直接返回
                if isinstance(resp, dict) and resp.get('code') == 'Disabled_Resource':
                    swanlog.warning(f"Experiment {http.exp_id} has been deleted, skipping column upload.")
                    return None
            raise e


__all__ = [
    "upload_logs",
    "upload_media_metrics",
    "upload_scalar_metrics",
    "upload_files",
    "upload_columns",
    "ScalarModel",
    "MediaModel",
    "ColumnModel",
    "FileModel",
]
