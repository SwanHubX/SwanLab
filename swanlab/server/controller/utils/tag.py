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
import os


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
    with open(file_path, "r") as f:
        lines = f.readlines()
        # 解析为list<dict>
        for i in range(len(lines)):
            lines[i] = json.loads(lines[i])
        return lines
