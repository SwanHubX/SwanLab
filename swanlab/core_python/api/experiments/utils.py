"""
@author: Zhou QiYang
@file: utils.py
@time: 2026/1/10 22:09
@description: 实验相关的后端API接口中的工具函数
"""

from swanlab.core_python.api.types import ColumnType


# 从前缀中获取指标类型
def parse_column_type(column: str) -> ColumnType:
    column_type = column.split('.', 1)[0]
    if column_type == 'summary':
        return 'SCALAR'
    elif column_type == 'config':
        return 'CONFIG'
    else:
        return 'STABLE'


# 将下划线命名转化为驼峰命名
def to_camel_case(name: str) -> str:
    return ''.join([w.capitalize() if i > 0 else w for i, w in enumerate(name.split('_'))])
