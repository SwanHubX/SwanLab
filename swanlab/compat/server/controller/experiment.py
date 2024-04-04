#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-02 21:23:30
@File: swanlab/compat/server/controller/experiment.py
@IDE: vscode
@Description:
    用于兼容实验部分的接口
"""
import shutil
from ....db import *
from urllib.parse import unquote
from swanlab.server.settings import (
    get_media_dir,
    get_tag_dir,
)
import os
import ujson


def _transfer_logs(log_path: str, text_path: str):
    """
    将原本存于某 txt 中的内容直接以字符串的形式存到对应的 log 文件
    一个 log 文件可能对应多个 txt 文件，该函数一次处理一个 log 文件

    Parameters
    ----------
    log_path : str
        tag 对应的 logs 下目录的路径
    text_path : str
        tag 对应 media 下目录的路径

    Returns
    -------
        单个文件转化后的数据
    """
    with open(log_path, "r") as f:
        data = f.read().split("\n")
        # 将每个元素转成dict
        data = [ujson.loads(line) for line in data if line]
    # 打开对应的 txt 文件
    for index, log in enumerate(data):
        path = os.path.join(text_path, log["data"])
        with open(path, "r") as f:
            text = f.read()
        data[index]["data"] = text
    # # 保存 log 文件
    with open(log_path, "w") as f:
        for line in data:
            f.write(ujson.dumps(line, ensure_ascii=False) + "\n")
    return data


def compat_text(experiment_id: int, tag: str):
    """
    0.3.0 之前，text 图表中的媒体数据单独存放于 txt 文件中
    但是因为其本质是字符串，可以直接存储于 logs 下

    Parameters
    ----------
    experiment_id : int
        实验 id
    tag : str
        tag 名

    Returns
    -------
        如果不需要转化，返回 False，如果需要转化，返回转化后的数据
    """

    tag = Tag.filter(Tag.experiment_id == experiment_id, Tag.name == unquote(tag)).first()
    if tag.type != "text":
        return False
    # 检测在 media 下是否存在同名目录
    run_id = Experiment.get(experiment_id).run_id
    text_path = get_media_dir(run_id, tag.name)
    # 是否存在该路径
    if not os.path.exists(text_path):
        return False
    # 将文件内容复制到 logs 记录中
    tag_path = get_tag_dir(run_id, tag.name)
    # 列出下面所有文件，除了 _summary.json
    files = [file for file in os.listdir(tag_path) if file != "_summary.json"]
    data: list = []
    for file in files:
        data.extend(_transfer_logs(os.path.join(tag_path, file), text_path))
    # 删除原来的媒体目录
    shutil.rmtree(text_path)
    return data
