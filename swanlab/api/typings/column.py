"""
@author: caddiesnew
@file: column.py
@time: 2026/4/20
@description: 公共查询 API 实验列类型定义
"""

from typing import Any, Dict, Optional, TypedDict

from .common import ApiColumnClassLiteral, ApiColumnDataTypeLiteral


class ApiColumnErrorType(TypedDict, total=False):
    """列错误信息"""

    message: str
    code: str


class ApiColumnType(TypedDict, total=False):
    """
    实验列数据类型

    注意：后端响应使用以下字段名：
    - class: 列的分类 (CUSTOM/SYSTEM)
    - type: 列的数据类型 (FLOAT/STRING/IMAGE等)
    - createdAt: 创建时间戳（蛇峰命名）
    """

    # 每个 column 与一个项目和实验绑定
    project_id: str
    run_id: str
    # 列的分类：CUSTOM 或 SYSTEM
    column_class: ApiColumnClassLiteral
    # 列的数据类型
    column_type: ApiColumnDataTypeLiteral
    # 列的键名（唯一标识）
    key: str
    # 列的显示名称，默认为 key 的值
    name: str
    # 创建时间戳
    createdAt: int
    # 错误信息
    error: Optional[Dict[str, Any]]


class ApiColumnCsvExportType(TypedDict):
    """列 CSV 导出响应类型"""

    # 临时下载 URL
    url: str
