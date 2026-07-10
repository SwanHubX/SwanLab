"""
@author: caddiesnew
@file: key.py
@time: 2026/7/9
@description: Key 实体类型定义 — House /metrics/*/keys 接口的请求/响应类型
"""

from typing import List, TypedDict


class ApiKeysPageType(TypedDict):
    """House /metrics/*/keys 单页响应。"""

    keys: List[str]
    nextCursor: str
    hasMore: bool
