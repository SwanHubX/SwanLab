"""
@author: Zhou Qiyang
@file: model.py
@time: 2025/12/18 20:10
@description: OpenApi 中的公用实体对象
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
        return str(self.name)
