"""
@author: Zhou QiYang
@file: projects.py
@time: 2026/1/10 21:48
@description: $END$
"""

from typing import TypedDict, List, Dict


# ------------------------------------- 通用类型 -------------------------------------
class ProjectLabelType(TypedDict):
    name: str  # 项目标签名称


# 项目信息
class ProjectType(TypedDict):
    cuid: str  # 项目CUID, 唯一标识符
    name: str  # 项目名
    path: str  # 项目路径
    url: str  # 项目URL
    description: str  # 项目描述
    visibility: str  # 可见性, 'PUBLIC' 或 'PRIVATE'
    createdAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    updatedAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    group: Dict[str, str]  # 包含项目所属工作空间名称 (workspace)
    projectLabels: List[ProjectLabelType]  # 项目标签
    _count: Dict[str, int]  # 项目的统计信息


# ------------------------------------- 后端返回信息 -------------------------------------
class ProjResponseType(TypedDict):
    list: List[ProjectType]  # 项目列表
    size: int  # 每页项目数量
    pages: int  # 总页数
    total: int  # 总项目数量
