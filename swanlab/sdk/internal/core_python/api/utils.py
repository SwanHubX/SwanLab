"""
@author: cunyue
@file: utils.py
@time: 2026/4/14 19:00
@description: API工具函数
"""

from swanlab.sdk.typings.run import SidebarItemType


def parse_column_type(column: str) -> SidebarItemType:
    """从前缀中获取指标类型"""
    column_type = column.split(".", 1)[0]
    if column_type == "summary":
        return "SCALAR"
    elif column_type == "config":
        return "CONFIG"
    else:
        return "STABLE"


def to_camel_case(name: str) -> str:
    """将下划线命名转化为驼峰命名"""
    return "".join([w.capitalize() if i > 0 else w for i, w in enumerate(name.split("_"))])
