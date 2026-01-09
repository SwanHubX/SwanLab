#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/8 14:36
@File: types.py
@IDE: pycharm
@Description:
    OpenAPI 相关数据结构
"""

from typing import Dict, Generic, List, TypeVar

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict


class BaseModel(PydanticBaseModel):
    def __getitem__(self, key):
        return getattr(self, key)


class Experiment(BaseModel):
    cuid: str               # 实验CUID, 唯一标识符
    name: str               # 实验名
    description: str = ""   # 实验描述
    state: str              # 实验状态, 'FINISHED' 或 'RUNNING'
    show: bool              # 显示状态
    createdAt: str          # e.g., '2024-11-23T12:28:04.286Z'
    finishedAt: str = ""    # e.g., '2024-11-23T12:28:04.286Z'
    user: Dict[str, str]    # 实验创建者, 包含 'username' 与 'name'
    profile: Dict           # 实验相关配置


class Project(BaseModel):
    cuid: str                   # 项目CUID, 唯一标识符
    name: str                   # 项目名
    description: str = ""       # 项目描述
    visibility: str             # 可见性, 'PUBLIC' 或 'PRIVATE'
    createdAt: str              # e.g., '2024-11-23T12:28:04.286Z'
    updatedAt: str              # e.g., '2024-11-23T12:28:04.286Z'
    group: Dict[str, str]       # 工作空间信息, 包含 'type', 'username', 'name'
    count: Dict[str, int] = {}  # 项目的统计信息


D = TypeVar("D")


class ApiResponse(BaseModel, Generic[D]):
    code: int    # HTTP状态码
    errmsg: str  # API错误消息, 只有请求错误时非空
    data: D      # 返回数据
    model_config = ConfigDict(arbitrary_types_allowed=True)


class Pagination(BaseModel, Generic[D]):
    total: int      # 总数
    list: List[D]   # 列表数据，泛型
