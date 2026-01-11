"""
@author: Zhou QiYang
@file: __init__.py.py
@time: 2026/1/11 23:44
@description: OpenApi 中的基础对象
"""

from typing import Dict

from .utils import Label


class ApiBase:
    @property
    def __dict__(self) -> Dict[str, object]:
        """
        Return a dictionary containing all @property fields.
        """
        result = {}
        cls = type(self)
        for attr_name in dir(cls):
            if attr_name.startswith('_'):
                continue
            attr = getattr(cls, attr_name, None)
            if isinstance(attr, property):
                result[attr_name] = self.__getattribute__(attr_name)
        return result


__all__ = ['ApiBase', 'Label']
