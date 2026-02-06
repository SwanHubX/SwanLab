"""
@author: Zhou Qiyang
@file: __init__.py
@time: 2026/1/27 20:43
@description: 工作空间
"""

from typing import Literal, Dict, Optional, List

from swanlab.api.utils import get_properties
from swanlab.core_python import Client
from swanlab.core_python.api.type import WorkspaceType, RoleType
from swanlab.core_python.auth.providers.api_key import LoginInfo


class Workspace:
    def __init__(
        self, client: Client, *, data: WorkspaceType, web_host: Optional[str] = None, login_info: Optional[LoginInfo] = None
    ) -> None:
        self._client = client
        self._data = data
        self._web_host = web_host
        self._login_info = login_info

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
        return self._data.get('profile', dict())

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

    def projects(
        self,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ):
        """
        Get all projects in this workspace.
        :param sort: Sort order, optional
        :param search: Search keyword, optional
        :param detail: Whether to return detailed info, optional
        :return: Projects instance, iterable to get project information
        """
        from swanlab.api.projects import Projects

        if self._web_host is None or self._login_info is None:
            raise RuntimeError("web_host and login_info are required. Use api.workspace() instead of creating Workspace directly.")
        return Projects(
            self._client,
            web_host=self._web_host,
            path=self._data['username'],
            sort=sort,
            search=search,
            detail=detail,
        )

    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)


__all__ = ['Workspace']
