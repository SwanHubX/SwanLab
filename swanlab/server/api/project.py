#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 00:38:23
@File: swanlab\server\api\project.py
@IDE: vscode
@Description:
    项目相关的api，前缀：/project
"""
from fastapi import APIRouter
from ..module.resp import SUCCESS_200
from ...database import PT

router = APIRouter()


# 列出当前项目下的所有实验
@router.get("")
async def _():
    """
    获取项目信息，列出当前项目下的所有实验
    """
    pt = PT()
    return SUCCESS_200(data=pt.get())
