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
from typing import List, Tuple
import asyncio
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
async def upload_logs(logs: List[str], level: str = "INFO"):
    """
    模拟一下，上传日志和实验指标信息
    :param logs: 日志列表
    :param level: 日志级别，'INFO', 'ERROR'，默认INFO
    """
    http = get_http()
    # 将logs解析为json格式
    metrics = [{"level": level, "message": x} for x in logs]
    data = create_data(metrics, "log")
    await http.post(url, data)
    await asyncio.sleep(1)


@async_error_handler
async def upload_media_metrics(media_metrics: List[Tuple[dict, str, str, List[str]]]):
    """
    上传指标的媒体数据
    :param media_metrics: 媒体指标数据，
        每个元素为元组，第一个元素为指标信息，
        第二个元素为指标的名称key，
        第三个元素为指标类型
        第四个元素为这个指标信息中包含的文件列表，每个元素为文件的绝对路径
    """
    print("上传媒体指标信息: ", media_metrics)


@async_error_handler
async def upload_scalar_metrics(scalar_metrics: List[dict]):
    """
    上传指标的标量数据
    """
    http = get_http()
    data = create_data(scalar_metrics, "scalar")
    await http.post(url, data)
    await asyncio.sleep(1)


@async_error_handler
async def upload_files(files: List[str]):
    """
    上传files文件夹中的内容
    :param files: 文件列表，内部为文件绝对路径
    """
    # 去重list
    files = list(set(files))
    files = [os.path.basename(x) for x in files]
    await asyncio.sleep(1)
    # print("上传文件信息: ", files)
