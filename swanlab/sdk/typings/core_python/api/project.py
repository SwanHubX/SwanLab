"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:02
@description: SwanLab 运行时项目API类型
"""

from typing import Dict, List, Literal, TypedDict


class _ProjectCount(TypedDict):
    # 项目历史实验数量
    experiments: int
    # 项目贡献者
    contributors: int
    # 项目协作者
    collaborators: int
    # 项目被clone次数
    clones: int


class ProjectLabelType(TypedDict):
    # 项目标签名称
    name: str


class ProjectType(TypedDict):
    # 项目名称
    name: str
    # 项目所属的用户名
    username: str
    # 项目路径 '/:username/:name'
    path: str
    # 项目可见性
    visibility: Literal["PUBLIC", "PRIVATE"]
    # 项目描述
    description: str
    # 项目所属工作空间
    group: Dict[str, str]
    # 项目标签
    projectLabels: List[ProjectLabelType]
    # 项目统计信息
    _count: _ProjectCount


class InitProjectType(TypedDict):
    # 项目名称
    name: str
    # 项目所属的用户名
    username: str
    # 项目路径 '/:username/:name'
    path: str


class ProjResponseType(TypedDict):
    # 项目列表
    list: List[ProjectType]
    # 每页项目数量
    size: int
    # 总页数
    pages: int
    # 总项目数量
    total: int
