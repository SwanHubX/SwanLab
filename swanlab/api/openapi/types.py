#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/8 14:36
@File: types.py
@IDE: pycharm
@Description:
    OpenAPI 相关数据结构
"""

from typing import TypedDict, Optional, Dict


class ExperimentProfile(TypedDict):
    config: Dict        # 实验配置参数
    metadata: Dict      # 实验元数据
    requirements: str   # 实验依赖项
    conda: str          # 实验的 Conda 环境信息

class Experiment(TypedDict):
    cuid: str   # 实验CUID, 唯一标识符
    name: str   # 实验名
    description: str    # 实验描述
    state: str          # 实验状态, 'FINISHED' 或 'RUNNING'
    show: str           # 显示状态, 'true' 或 'false'
    createdAt: str              # e.g., '2024-11-23T12:28:04.286Z'
    finishedAt: Optional[str]   # e.g., '2024-11-23T12:28:04.286Z', 若不存在则为 None
    user: Dict[str, str]        # 实验创建者, 包含 'username' 与 'name'
    profile: ExperimentProfile  # 实验相关配置

class ApiErrorResponse(TypedDict):
    code: int       # HTTP错误代码
    message: str    # API错误消息
