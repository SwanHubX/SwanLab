#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 16:06:49
@File: swanlab/data/run/utils.py
@IDE: vscode
@Description:
    运行时工具函数
"""
from ...utils.file import check_tag_format, get_a_lock, check_exp_name_format, check_desc_format
from ...utils import get_package_version, create_time, generate_color
import datetime


def json_serializable(obj: dict):
    """
    Convert an object into a JSON-serializable format.
    """
    # 如果对象是基本类型，则直接返回
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj

    # 将日期和时间转换为字符串
    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()

    # 对于列表和元组，递归调用此函数
    if isinstance(obj, (list, tuple)):
        return [json_serializable(item) for item in obj]

    # 对于字典，递归调用此函数处理值
    if isinstance(obj, dict):
        return {key: json_serializable(value) for key, value in obj.items()}

    # 对于其他不可序列化的类型，转换为字符串表示
    return str(obj)
