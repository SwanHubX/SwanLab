"""
@author: Zhou QiYang
@file: project.py
@time: 2026/1/5 17:58
@description: OpenApi 中的项目对象
"""

from typing import List, Dict

from swanlab.api.utils import ApiBase, Label
from swanlab.core_python.api.type import ProjectType


class Project(ApiBase):
    """
    Representing a single project with some of its properties.
    """

    def __init__(self, data: ProjectType, web_host: str) -> None:
        self._data = data
        self._web_host = web_host

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

    @property
    def workspace(self) -> str:
        """
        Project workspace name.
        """
        return self._data["group"]["username"]

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
