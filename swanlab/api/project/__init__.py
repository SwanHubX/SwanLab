"""
@author: Zhou QiYang
@file: project.py
@time: 2026/1/5 17:58
@description: OpenApi 中的项目对象
"""

from functools import cached_property
from typing import List, Dict, Optional

from swanlab.api.experiments import Experiments
from swanlab.api.utils import Label, get_properties
from swanlab.api.workspace import Workspace
from swanlab.core_python.api.type import ProjectType
from swanlab.core_python.api.user import get_workspace_info
from swanlab.core_python.auth.providers.api_key import LoginInfo
from swanlab.core_python.client import Client


class Project:
    """
    Representing a single project with some of its properties.
    """

    def __init__(
        self, client: Client, *, web_host: str, data: ProjectType, login_info: Optional[LoginInfo] = None
    ) -> None:
        self._client = client
        self._web_host = web_host
        self._data = data
        self._login_info = login_info

    @property
    def name(self) -> str:
        """
        Project name.
        """
        return self._data['name']

    @property
    def path(self) -> str:
        """
        Project path in the format 'username/project-name'.
        """
        return self._data['path']

    @property
    def url(self) -> str:
        """
        Full URL to access the project.
        """
        return f"{self._web_host}/@{self._data['path']}"

    @property
    def description(self) -> str:
        """
        Project description.
        """
        return self._data['description']

    @property
    def visibility(self) -> str:
        """
        Project visibility, either 'PUBLIC' or 'PRIVATE'.
        """
        return self._data['visibility']

    @property
    def created_at(self) -> str:
        """
        Project creation timestamp
        """
        return self._data['createdAt']

    @property
    def updated_at(self) -> str:
        """
        Project last update timestamp
        """
        return self._data['updatedAt']

    @cached_property
    def workspace(self) -> Workspace:
        """
        Project workspace object.
        """
        data = get_workspace_info(self._client, path=self._data["group"]["username"])
        return Workspace(self._client, data=data, web_host=self._web_host, login_info=self._login_info)

    @property
    def labels(self) -> List[Label]:
        """
        List of Label attached to this project.
        """
        return [Label(label['name']) for label in self._data['projectLabels']]

    @property
    def count(self) -> Dict[str, int]:
        """
        Project statistics dictionary containing:
        experiments, contributors, children, collaborators, runningExps.
        """
        return self._data['_count']

    def runs(self, filters: Optional[Dict[str, object]] = None) -> Experiments:
        """
        Get all experiments in this project.
        :param filters: Filter conditions for experiments, optional
        :return: Experiments instance, iterable to get experiment information
        """
        if self._login_info is None:
            raise RuntimeError("login_info is required to access runs. Use api.project() instead of creating Project directly.")
        return Experiments(self._client, path=self._data['path'], login_info=self._login_info, filters=filters)

    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)