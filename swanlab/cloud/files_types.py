#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/5 16:55
@File: files_types.py
@IDE: pycharm
@Description:
    文件资源类型，分为四种：
    1. 日志类型，字符串方式解析
    2. 指标类型，JSON方式解析
    3. 文件类型，对应yml、json、txt等格式进行解析
    4. 媒体类型，文件类型解析
在本模块针对上述四种方式定义不同的类和不同的处理方式，每种类型对应一个数据上传接口
"""
from enum import Enum
from typing import List
from swanlab.error import NetworkError
from requests.exceptions import RequestException
from swanlab.log import swanlog


def async_error_handler(func):
    """
    用于进行统一的错误捕获
    """

    async def wrapper(*args, **kwargs):
        try:
            # 在装饰器中调用被装饰的异步函数
            result = await func(*args, **kwargs)
            return result
        except RequestException:
            swanlog.error('network error, swanlab will resume uploads when the network improves')
            return NetworkError()
        except Exception as e:
            return e

    return wrapper


@async_error_handler
async def mock_upload_logs(logs: List[str]):
    """
    模拟一下，上传日志和实验指标信息
    """
    print("上传日志信息: ", logs)


@async_error_handler
async def mock_upload_metrics(metrics: List[dict]):
    """
    模拟一下，上传日志和实验指标信息
    """
    print("上传指标信息: ", metrics)


@async_error_handler
async def mock_upload_files(files: List[str]):
    """
    模拟一下，上传日志和实验指标信息
    :param files: 文件列表，内部为文件绝对路径
    """
    print("上传文件信息: ", files)
    raise RequestException()


@async_error_handler
async def mock_upload_media(media: List[str]):
    """
    模拟一下，上传日志和实验指标信息
    :return:
    """
    print("上传媒体信息: ", media)


class FileType(Enum):
    """
    文件类型枚举，在此处定义文件类型以及不同的处理方式
    priority属性表示优先级，数字越大优先级越高，代表在一次数据上传的循环中先上传
    {
        "name": 枚举类型名称,
        "upload": 上传数据的方法,
        "priority": 优先级,
        "error": 错误描述
    }
    """
    LOG = {
        "name": "log",
        "upload": mock_upload_logs,
        "priority": 1,
        "error": "logs upload failed"
    }
    METRIC = {
        "name": "metric",
        "upload": mock_upload_metrics,
        "priority": 2,
        "error": "metrics upload failed"
    }
    FILE = {
        "name": "file",
        "upload": mock_upload_files,
        "priority": 1,
        "error": "files upload failed"
    }
    MEDIA = {
        "name": "media",
        "upload": mock_upload_media,
        "priority": 9,  # 最高优先级
        "error": "media upload failed"
    }
