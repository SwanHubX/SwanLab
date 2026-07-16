"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:02
@description: SwanLab 运行时项目API类型
"""

from typing import Literal, TypedDict

from typing_extensions import NotRequired


class _ProjectCount(TypedDict):
    # 项目历史实验数量
    experiments: int
    # 项目贡献者
    contributors: int
    # 项目协作者
    collaborators: int
    # 项目被clone次数
    clones: int


class _ProjectGroup(TypedDict):
    username: str


class ProjectType(TypedDict):
    # 项目ID
    cuid: str
    # 项目名称
    name: str
    # 项目版本，如果不存在此字段则按照最低版本处理
    version: NotRequired[Literal[1]]
    # 项目路径 '/:username/:name'
    path: str
    # 项目可见性
    visibility: Literal["PUBLIC", "PRIVATE"]
    # 项目所属空间
    group: _ProjectGroup
    # 项目统计信息
    _count: _ProjectCount
