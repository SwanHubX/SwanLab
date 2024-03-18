#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-17 15:13:31
@File: swanlab/server/controller/utils/tag.py
@IDE: vscode
@Description:
    tag相关处理函数
"""
from typing import List
import json
import ujson
import os

# tag 总结文件名
TAG_SUMMARY_FILE = "_summary.json"
# logs 目录下的配置文件
LOGS_CONFIGS = [TAG_SUMMARY_FILE]


def read_tag_data(file_path: str) -> List[dict]:
    """读取某一个tag文件数据，数据依据其中的index字段排序

    Parameters
    ----------
    file_path : str
        tag文件路径，绝对路径
    """
    # 如果文件内容为空，返回空列表
    if os.path.getsize(file_path) == 0:
        return []
    # 读取文件内容，文件内部本质上是字符串一堆json格式的数据，每一行是一个json并且有换行符分隔，需要拿到然后解析为list<dict>
    data = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        # 解析为list<dict>
        for i in range(len(lines)):
            if len(lines[i]):
                data.append(json.loads(lines[i]))
        return data


def get_tag_files(tag_path: str, exclude: List[str] = []) -> List[str]:
    """
    获取实验数据，并且做向下兼容，在v0.2.4版本以前的实验日志数据将转换为新的格式

    Parameters
    ----------
    file_path : str
        tag的路径
    exclude : List[str]
        需要排除的文件列表

    Returns
    -------
    List[str]
        tag文件列表
    """
    # 降序排列，最新的数据在最前面
    files: list = os.listdir(tag_path)
    previous_logs = [f for f in files if f.endswith(".json") and f not in exclude]
    current_logs = [f for f in files if f.endswith(".log")]
    current_logs.sort()
    # COMPAT 如果目标文件夹不存在*.log文件但存在*.json文件，说明是之前的日志格式，需要转换
    if len(current_logs) == 0 and len(previous_logs) > 0:
        for file in previous_logs:
            with open(os.path.join(tag_path, file), "r") as f:
                data = ujson.load(f)
                with open(os.path.join(tag_path, file.replace(".json", ".log")), "a") as f:
                    for d in data["data"]:
                        f.write(ujson.dumps(d) + "\n")
        current_logs = [f for f in os.listdir(tag_path) if f.endswith(".log")]
        current_logs.sort()
    return current_logs
