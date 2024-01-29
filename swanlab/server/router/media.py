#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-29 14:14:02
@File: swanlab\server\router\media.py
@IDE: vscode
@Description:
    多媒体相关接口路由
"""


from fastapi import APIRouter

from ..controller.media import (
    # 解析、获取音频数据
    get_audio_data,
)

router = APIRouter()


# ---------------------------------- 音频相关 ----------------------------------


@router.get("/audio")
def _(path: str):
    return get_audio_data(path)
