"""
@author: caddiesnew
@file: common.py
@time: 2026/4/20
@description: 公共查询 API 通用类型定义
"""

from typing import List, TypedDict


class ApiLabelType(TypedDict):
    name: str


class ApiPaginationType(TypedDict):
    list: List
    size: int
    pages: int
    total: int
