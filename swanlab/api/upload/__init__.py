#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:56
@File: __init__.py
@IDE: pycharm
@Description:
    上传相关接口
"""
from ..http import get_http, async_error_handler
from typing import List, Tuple, Dict
import json
import yaml
import os

url = '/house/metrics'


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


@async_error_handler
async def upload_logs(logs: List[dict], level: str = "INFO"):
    """
    模拟一下，上传日志和实验指标信息
    :param logs: 日志列表
    :param level: 日志级别，'INFO', 'ERROR'，默认INFO
    """
    http = get_http()
    # 将logs解析为json格式
    metrics = [{"level": level, **x} for x in logs]
    data = create_data(metrics, "log")
    await http.post(url, data)


@async_error_handler
async def upload_media_metrics(media_metrics: List[Tuple[dict, str, str, str]]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据，
        每个元素为元组，第一个元素为指标信息，
        第二个元素为指标的名称key，经过URI编码
        第三个元素为指标类型
        第四个元素为media文件夹路径
    """
    http = get_http()
    # 需要上传的文件路径[key, local_path]
    file_paths: Dict[str, str] = {}
    for metric, key, data_type, media_folder in media_metrics:
        if data_type == "text":
            # 字符串类型没有文件路径
            continue
        if isinstance(metric["data"], str):
            local_path = metric["data"]
            metric["data"] = "{}/{}".format(key, metric["data"])
            # 将文件路径添加到files_path中
            file_paths[metric["data"]] = os.path.join(media_folder, key, local_path)
        else:
            local_paths = metric['data']
            metric['data'] = ["{}/{}".format(key, x) for x in local_paths]
            for i, local_path in enumerate(local_paths):
                file_paths[metric['data'][i]] = os.path.join(media_folder, key, local_path)
    # 上传文件，先上传资源文件，再上传指标信息
    keys = list(file_paths.keys())
    local_paths = list(file_paths.values())
    result = await http.upload_files(keys, local_paths)
    # result = await http.upload(keys[0], local_paths[0])
    # 上传指标信息
    await http.post(url, create_data([x[0] for x in media_metrics], "media"))


@async_error_handler
async def upload_scalar_metrics(scalar_metrics: List[dict]):
    """
    上传指标的标量数据
    """
    http = get_http()
    data = create_data(scalar_metrics, "scalar")
    await http.post(url, data)


_valid_files = {
    'config.yaml': ['config', 'yaml'],
    'requirements.txt': ['requirements', 'txt'],
    'swanlab-metadata.json': ['metadata', 'json']
}
"""
支持上传的文件列表，filename: key
"""


@async_error_handler
async def upload_files(files: List[str]):
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
        with open(filepath, 'r') as f:
            if _valid_files[filename][1] == 'json':
                data[_valid_files[filename][0]] = json.load(f)
            elif _valid_files[filename][1] == 'yaml':
                data[_valid_files[filename][0]] = yaml.load(f, Loader=yaml.FullLoader)
            else:
                data[_valid_files[filename][0]] = f.read()

    await http.put(f'/project/{http.groupname}/{http.projname}/runs/{http.exp_id}/profile', data)
