"""
@author: Zhou QiYang
@file: _user.py
@time: 2026/1/11 16:57
@description: $END$
"""

from swanlab.api.base import ApiBase
from swanlab.core_python.api.experiments.type import UserType


class User(ApiBase):
    def __init__(self, data: UserType) -> None:
        self._data = data

    @property
    def name(self) -> str:
        return self._data['name']

    @property
    def username(self) -> str:
        return self._data['username']
