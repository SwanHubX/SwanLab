"""
@author: Zhou QiYang
@file: __init__.py
@time: 2026/1/11 23:44
@description: OpenApi 中的基础对象
"""

from dataclasses import dataclass


@dataclass
class Label:
    """
    Project label object
    you can get the label name by str(label)
    """

    name: str

    def __str__(self) -> str:
        return self.name


__all__ = ['Label']
