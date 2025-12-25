"""
@author: Zhou Qiyang
@file: types.py
@time: 2025/12/17 16:35
@description: OpenApi 用到的类型文件
"""

from typing import TypedDict, Optional, List, Dict, Literal


# ------------------------------------- 通用类型 -------------------------------------
class GroupType(TypedDict):
    username: str  # 工作空间名称 (workspace)


class ProjectLabelType(TypedDict):
    name: str  # 项目标签名称


class UserType(TypedDict):
    username: str
    name: str


StateType = Literal['FINISHED', 'CRASHED', 'ABORTED', 'RUNNING']


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
    group: GroupType  # 项目所属工作空间名称 (workspace)
    projectLabels: List[ProjectLabelType]  # 项目标签
    _count: Dict[str, int]  # 项目的统计信息


class RunType(TypedDict):
    cuid: str  # 实验CUID, 唯一标识符
    name: str  # 实验名称
    createdAt: str  # 创建时间, e.g., '2024-11-23T12:28:04.286Z'
    description: str  # 实验描述
    labels: List[ProjectLabelType]  # 实验标签列表
    profile: Dict[str, Dict[str, object]]  # 实验配置和摘要信息，包含 'config' 和 'scalar'
    state: StateType  # 实验状态
    cluster: str  # 实验组
    job: str  # 任务类型
    runtime: str  # 运行时间
    user: UserType  # 实验所属用户


# ------------------------------------- 发送到后端请求体中的参数 -------------------------------------
class ProjParamType(TypedDict):
    page: int  # 页码
    size: int  # 每页项目数量
    sort: Optional[List[str]]  # 排序方式（包含多个条件的列表）
    search: Optional[str]  # 搜索关键词
    detail: Optional[bool]  # 是否返回详细信息（_count）


# ------------------------------------- 后端返回信息 -------------------------------------
class ProjResponseType(TypedDict):
    list: List[ProjectType]  # 项目列表
    size: int  # 每页项目数量
    pages: int  # 总页数
    total: int  # 总项目数量
