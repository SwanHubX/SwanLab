"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:02
@description: SwanLab 运行时项目API类型
"""

from typing import TypedDict


class InitProjectType(TypedDict):
    # 项目名称
    name: str
    # 项目所属的用户名
    username: str
    # 项目路径 '/:username/:name'
    path: str
