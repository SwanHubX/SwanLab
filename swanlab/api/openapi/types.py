#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/8 14:36
@File: types.py
@IDE: pycharm
@Description:
    OpenAPI 相关数据结构
"""

from typing import Optional, Dict, TypeVar, Generic, List, Any
from pydantic import BaseModel, ConfigDict


class Experiment(BaseModel):
    cuid: str  # 实验CUID, 唯一标识符
    name: str  # 实验名
    description: Optional[str]  # 实验描述
    state: str  # 实验状态, 'FINISHED' 或 'RUNNING'
    show: bool  # 显示状态
    createdAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    finishedAt: Optional[str]  # e.g., '2024-11-23T12:28:04.286Z', 若不存在则为 None
    user: Dict[str, str]  # 实验创建者, 包含 'username' 与 'name'
    profile: Dict  # 实验相关配置


class Project(BaseModel):
    cuid: str  # 项目CUID, 唯一标识符
    name: str  # 项目名
    description: Optional[str]  # 项目描述
    visbility: str  # 可见性, 'PUBLIC' 或 'PRIVATE'
    createdAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    updatedAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    path: str  # 项目路径, e.g., '/project/username/project_name'
    group: Dict[str, Any]  # 工作空间信息, 包含 'type', 'username', 'name'

    _count: Optional[Dict[str, Any]]  # 仅当 detail=True 时返回, 包含项目的详细信息

    model_config = ConfigDict(extra='allow', validate_assignment=True)  # 显式处理 _count 字段


D = TypeVar("D")


class ApiResponse(BaseModel, Generic[D]):
    code: int  # HTTP状态码
    errmsg: str  # API错误消息, 只有请求错误时非空
    data: D  # 返回数据


class Pagination(BaseModel, Generic[D]):
    total: int  # 总数
    list: List[D]  # 列表数据，类型为泛型 T
