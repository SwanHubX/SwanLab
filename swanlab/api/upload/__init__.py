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


@async_error_handler
async def upload_logs(logs: List[str], level: str = "INFO"):
    """
    模拟一下，上传日志和实验指标信息
    :param logs: 日志列表
    :param level: 日志级别，'INFO', 'ERROR'，默认INFO
    """
    print("上传日志信息: ", logs)


@async_error_handler
async def upload_media_metrics(media_metrics: List[Tuple[dict, List[str]]]):
    """
    上传指标的媒体数据
    """
    print("上传媒体指标信息: ", media_metrics)


@async_error_handler
async def upload_scalar_metrics(scalar_metrics: List[dict]):
    """
    上传指标的标量数据
    """
    print("上传标量指标信息: ", scalar_metrics)


@async_error_handler
async def upload_files(files: List[str]):
    """
    模拟一下，上传日志和实验指标信息
    :param files: 文件列表，内部为文件绝对路径
    """
    # 去重list
    files = list(set(files))
    files = [os.path.basename(x) for x in files]
    await asyncio.sleep(10)
    # print("上传文件信息: ", files)
