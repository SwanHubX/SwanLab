#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:56
@File: __init__.py
@IDE: pycharm
@Description:
    上传相关接口
"""
from ..http import get_http, sync_error_handler
from .model import ColumnModel, MediaModel, ScalarModel
from typing import List, Dict, ByteString
from swanlab.error import FileError, ApiError
from swanlab.log import swanlog
import json
import yaml
import os

house_url = '/house/metrics'


def create_data(metrics: List[dict], metrics_type: str) -> dict:
    """
    携带上传日志的指标信息
    """
    http = get_http()
    return {
        "projectId": http.proj_id,
        "experimentId": http.exp_id,
        "type": metrics_type,
        "metrics": metrics
    }


@sync_error_handler
def upload_logs(logs: List[dict], level: str = "INFO"):
    """
    上传日志信息
    :param logs: 日志列表
    :param level: 日志级别，'INFO', 'ERROR'，默认INFO
    """
    http = get_http()
    # 将logs解析为json格式
    metrics = [{"level": level, **x} for x in logs]
    data = create_data(metrics, "log")
    http.post(house_url, data)


@sync_error_handler
def upload_media_metrics(media_metrics: List[MediaModel]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据集合
    """
    http = get_http()
    # 需要上传的文件[path, raw]
    file_paths: Dict[str, ByteString] = {}
    metrics = []
    for media in media_metrics:
        metric = media.to_dict()
        metrics.append(metric)
        if media.raw is None:
            continue
        key_encoded = media.key_encoded
        file_names: List[str] = metric["data"]
        for file_name in file_names:
            file_paths["{}/{}".format(key_encoded, file_name)] = media.raw
    # 上传文件，先上传资源文件，再上传指标信息
    keys = list(file_paths.keys())
    local_paths = list(file_paths.values())
    http.upload_files(keys, local_paths)
    # 上传指标信息
    http.post(house_url, create_data([x[0] for x in media_metrics], "media"))


@sync_error_handler
def upload_scalar_metrics(scalar_metrics: List[ScalarModel]):
    """
    上传指标的标量数据
    """
    http = get_http()
    data = create_data([x.to_dict() for x in scalar_metrics], ScalarModel.type)
    http.post(house_url, data)


_valid_files = {
    'config.yaml': ['config', 'yaml'],
    'requirements.txt': ['requirements', 'txt'],
    'swanlab-metadata.json': ['metadata', 'json']
}
"""
支持上传的文件列表，filename: key
"""


@sync_error_handler
def upload_files(files: List[str]):
    """
    上传files文件夹中的内容
    :param files: 文件列表，内部为文件绝对路径
    """
    http = get_http()
    # 去重list
    files = list(set(files))
    files = {os.path.basename(x): x for x in files}
    # 读取文件配置，生成更新信息
    data = {}
    for filename, filepath in files.items():
        if filename not in _valid_files:
            continue
        try:
            with open(filepath, 'r') as f:
                if _valid_files[filename][1] == 'json':
                    data[_valid_files[filename][0]] = json.load(f)
                elif _valid_files[filename][1] == 'yaml':
                    d = yaml.load(f, Loader=yaml.FullLoader)
                    if d is None:
                        raise FileError
                    data[_valid_files[filename][0]] = d
                else:
                    data[_valid_files[filename][0]] = f.read()
        except json.decoder.JSONDecodeError:
            raise FileError
    http.put(f'/project/{http.groupname}/{http.projname}/runs/{http.exp_id}/profile', data)


@sync_error_handler
def upload_column(columns: List[ColumnModel]):
    """
    上传列信息，需要注意的是一次只能上传一个列，所以函数名不带s
    但是在设计上是聚合上传的，所以在此处需要进行拆分然后分别上传
    """
    http = get_http()
    url = f'/experiment/{http.exp_id}/column'
    # WARNING 这里不能使用并发请求，可见 https://github.com/SwanHubX/SwanLab-Server/issues/113
    for column in columns:
        try:
            http.post(url, column.to_dict())
        except ApiError as e:
            swanlog.error(f"Upload column {column.key} failed: {e.resp.status_code}")


__all__ = [
    "upload_logs",
    "upload_media_metrics",
    "upload_scalar_metrics",
    "upload_files",
    "upload_column",
    "ScalarModel",
    "MediaModel",
    "ColumnModel"
]
