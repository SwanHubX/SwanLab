#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-29 14:14:02
@File: swanlab\server\router\static.py
@IDE: vscode
@Description:
    多媒体相关接口路由
    暂时只有一个，用于获取媒体数据，这媒体包括音频、视频、图片、文字等，返回的接口是统一的
    返回的数据格式也不再是JSON，而是二进制数据，这样可以减少数据传输的大小
    我们约定前端通过一个相对路径获取媒体数据，这个路径相对于当前watch的文件夹
"""

from ..settings import get_media_dir
from fastapi import APIRouter
from urllib.parse import quote
from fastapi.responses import FileResponse
from ...db import Experiment
import os

router = APIRouter()


# ---------------------------------- 音频相关 ----------------------------------


@router.get("/{path:path}")
def _(path: str, tag: str, experiment_id: str):
    """获取媒体文件
    通过参数拼接为本地路径，返回二进制数据
    tag需要进行url编码
    """
    run_id = Experiment.get(Experiment.id == experiment_id).run_id
    media_path = os.path.join(get_media_dir(run_id, quote(tag, safe="")), path)
    return FileResponse(media_path)
