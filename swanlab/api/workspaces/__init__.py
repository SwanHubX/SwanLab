"""
@author: Zhou Qiyang
@file: __init__.py
@time: 2026/1/27 18:13
@description: 工作空间迭代器
"""

from typing import Iterator

from swanlab.api.workspace import Workspace
from swanlab.core_python import Client
from swanlab.core_python.api.user import get_user_groups, get_workspace_info


class Workspaces:
    def __init__(self, client: Client, *, username: str) -> None:
        self._client = client
        self._username = username

    def get_all_workspaces(self, username: str = None):
        """Get all workspaces of specific user (defaults to current user)"""
        cur_username = username if username else self._username
        resp = get_user_groups(self._client, username=cur_username)
        groups = [r['username'] for r in resp]
        return [cur_username] + groups

    def __iter__(self) -> Iterator[Workspace]:
        for space in self.get_all_workspaces():
            data = get_workspace_info(self._client, workspace=space)
            yield Workspace(data=data)
