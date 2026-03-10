"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:02
@description: SwanLab 运行时项目API类型
"""

from typing import Literal, TypedDict


class _ProjectCount(TypedDict):
    # 项目历史实验数量
    experiments: int
    # 项目贡献者
    contributors: int
    # 项目协作者
    collaborators: int
    # 项目被clone次数
    clones: int


class ProjectType(TypedDict):
    # 项目名称
    name: str
    # 项目所属的用户名
    username: str
    # 项目路径 '/:username/:name'
    path: str
    # 项目可见性
    visibility: Literal["PUBLIC", "PRIVATE"]
    # 项目统计信息
    _count: _ProjectCount


class InitProjectType(TypedDict):
    # 项目名称
    name: str
    # 项目所属的用户名
    username: str
    # 项目路径 '/:username/:name'
    path: str
