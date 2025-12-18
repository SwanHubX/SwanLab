"""
@author: Zhou Qiyang
@file: types.py
@time: 2025/12/17 16:35
@description: OpenApi 用到的类型文件
"""
from typing import TypedDict, Optional, List, Dict


# 发送到后端查询项目信息的字段
class ProjParamType(TypedDict):
    page: int  # 页码
    size: int  # 每页项目数量
    sort: Optional[List[str]]  # 排序方式（包含多个条件的列表）
    search: Optional[str]  # 搜索关键词
    detail: Optional[bool]  # 是否返回详细信息（_count）


class GroupType(TypedDict):
    username: str  # 组织或者个人的名称标识
    status: str  # 组织或者个人的状态
    type: str  # 用户身份（个人、团队）


class ProjectLabelType(TypedDict):
    name: str  # 项目标签名称


class ProjectType(TypedDict):
    cuid: str  # 项目CUID, 唯一标识符
    name: str  # 项目名
    path: str  # 项目路径
    url: str  # 项目URL
    description: str  # 项目描述
    visibility: str  # 可见性, 'PUBLIC' 或 'PRIVATE'
    createdAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    updatedAt: str  # e.g., '2024-11-23T12:28:04.286Z'
    projectLabels: List[ProjectLabelType]  # 项目标签
    group: GroupType  # 团队信息
    _count: Dict[str, int]  # 项目的统计信息


# 后端返回的项目信息
class ProjResponseType(TypedDict):
    list: List[ProjectType]  # 项目列表
    size: int  # 每页项目数量
    pages: int  # 总页数
    total: int  # 总项目数量
