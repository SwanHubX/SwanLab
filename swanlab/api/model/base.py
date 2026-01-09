"""
@author: Zhou Qiyang
@file: model.py
@time: 2025/12/18 20:10
@description: OpenApi 中的基础对象
"""

from typing import Dict

from swanlab.core_python.api.type import ProjectLabelType, UserType


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


class Label(ApiBase):
    """
    Project label object
    you can get the label name by str(label)
    """

    def __init__(self, data: ProjectLabelType) -> None:
        self._data = data

    @property
    def name(self) -> str:
        """
        Label name.
        """
        return self._data['name']

    def __str__(self) -> str:
        return str(self.name)


class User(ApiBase):
    def __init__(self, data: UserType) -> None:
        self._data = data

    @property
    def name(self) -> str:
        return self._data['name']

    @property
    def username(self) -> str:
        return self._data['username']
