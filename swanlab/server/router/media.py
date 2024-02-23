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
from typing import List
import asyncio
import os
from ..module.resp import SUCCESS_200

router = APIRouter()


# ---------------------------------- 音频相关 ----------------------------------


@router.get("/audio/{path:path}")
def _(path: str, tag: str, run_id: str):
    """获取媒体文件
    通过参数拼接为本地路径，返回二进制数据
    tag需要进行url编码
    """
    media_path = os.path.join(get_media_dir(run_id, quote(tag, safe="")), path)
    return FileResponse(media_path)


# ---------------------------------- 文字相关 ----------------------------------


def read_file_sync(file_name):
    with open(file_name, "r", encoding="utf8") as file:
        return file.read()


async def read_file_async(file_name):
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(None, read_file_sync, file_name)


@router.get("/text")
async def _(path, tag: str, run_id: str):
    """获取文本内容"""

    paths = path.split(",") if "," in path else [path]

    tasks = [read_file_async(os.path.join(get_media_dir(run_id, quote(tag, safe="")), p)) for p in paths]
    res = await asyncio.gather(*tasks)

    return SUCCESS_200({"text": res})
