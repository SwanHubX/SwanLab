"""
@author: Zhou Qiyang
@file: __init__.py
@time: 2026/1/27 20:43
@description: 工作空间
"""

from typing import Literal, Dict

from swanlab.api.utils import get_properties
from swanlab.core_python import Client
from swanlab.core_python.api.type import WorkspaceType, RoleType


class Workspace:
    def __init__(self, client: Client, *, data: WorkspaceType):
        self._client = client
        self._data = data

    @property
    def name(self) -> str:
        """
        Workspace display name.
        """
        return self._data['name']

    @property
    def username(self) -> str:
        """
        Workspace name.
        """
        return self._data['username']

    @property
    def workspace_type(self) -> Literal['TEAM', 'PERSON']:
        """
        Workspace type.
        """
        return self._data['type']

    @property
    def profile(self) -> Dict[str, str]:
        """
        Workspace profile.
        """
        return self._data['profile']

    @property
    def comment(self) -> str:
        """
        Workspace comment.
        """
        return self._data['comment']

    @property
    def role(self) -> RoleType:
        """
        Current login user's role in the workspace (only display when type=TEAM).
        """
        return self._data['role']

    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)


__all__ = ['Workspace']
