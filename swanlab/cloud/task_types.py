#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/5 16:55
@File: task_types.py
@IDE: pycharm
@Description:
    文件资源类型，分为四种：
    1. 日志类型，字符串方式解析
    2. 标量指标类型，JSON方式解析，上传指标信息
    3. 媒体指标类型，JSON+文件方式解析，同时上传资源文件和指标信息
    3. 文件类型，对应yml、json、txt等格式进行解析
在本模块针对上述四种方式定义不同的类和不同的处理方式，每种类型对应一个数据上传接口
"""
from enum import Enum
from swanlab.api.upload import *


class UploadType(Enum):
    """
    上传类型枚举，在此处定义不同的处理方式
    """
    LOG = {
        "upload": upload_logs,
    }
    """
    上传日志字符串
    """

    SCALAR_METRIC = {
        "upload": upload_scalar_metrics,
    }
    """
    上传标量指标
    """

    MEDIA_METRIC = {
        "upload": upload_media_metrics,
    }
    """
    上传媒体指标
    """

    FILE = {
        "upload": upload_files,
    }
    """
    上传实验信息配置文件
    """

    COLUMN = {
        # 一次只能上传一个列，所以函数名不带s
        "upload": upload_column,
    }
    """
    上传列信息
    """
