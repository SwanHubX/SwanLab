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
from ..utils import ResponseBody
from ...env import SWANLAB_LOGS_FOLDER
import os
import ujson

router = APIRouter()


# 列出当前项目下的所有实验
@router.get("/experiments")
async def _():
    config_path = os.path.join(SWANLAB_LOGS_FOLDER, "project.json")
    experiments = ujson.load(open(config_path, "r"))
    return ResponseBody(0, data=experiments)
