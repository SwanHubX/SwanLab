#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/8 14:36
@File: types.py
@IDE: pycharm
@Description:
    OpenAPI 相关数据结构
"""

from typing import Optional, Dict, TypeVar, Generic, List
from pydantic import BaseModel

class Experiment(BaseModel):
    cuid: str   # 实验CUID, 唯一标识符
    name: str   # 实验名
    description: Optional[str]    # 实验描述
    state: str          # 实验状态, 'FINISHED' 或 'RUNNING'
    show: bool           # 显示状态
    createdAt: str              # e.g., '2024-11-23T12:28:04.286Z'
    finishedAt: Optional[str]   # e.g., '2024-11-23T12:28:04.286Z', 若不存在则为 None
    user: Dict[str, str]        # 实验创建者, 包含 'username' 与 'name'
    profile: Dict  # 实验相关配置

D = TypeVar("D")

class ApiResponse(BaseModel, Generic[D]):
    code: int       # HTTP状态码
    errmsg: str    # API错误消息, 只有请求错误时非空
    data: D  # 返回数据

class Pagination(BaseModel, Generic[D]):
    total: int  # 总数
    list: List[D]  # 列表数据，类型为泛型 T
