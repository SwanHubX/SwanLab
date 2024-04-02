#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-02 21:23:30
@File: swanlab/compat/server/controller/experiment.py
@IDE: vscode
@Description:
    用于兼容实验部分的接口
"""
from ....db import *
from urllib.parse import unquote
from swanlab.server.settings import (
    get_media_dir,
    get_tag_dir,
)
import os


def compat_text(experiment_id: int, tag: str):
    """
    0.3.0 之前，text 图表中的媒体数据单独存放于 txt 文件中
    但是因为其本质是字符串，可以直接存储于 logs 下
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
    print(files)
    return True
