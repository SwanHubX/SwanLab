"""
@author: Zhou Qiyang
@file: model.py
@time: 2025/12/18 20:10
@description: OpenApi查询结果将以对象返回，并且对后端的返回字段进行一些筛选
"""

from typing import List, Dict, Optional

from swanlab.core_python import Client
from swanlab.core_python.api.project import get_entity_projects
from .type import ProjectType, ProjectLabelType


class Label:
    """
    Project label object
    you can get the label name by str(label)
    """

    def __init__(self, data: ProjectLabelType):
        self._data = data

    @property
    def name(self):
        """
        Label name.
        """
        return self._data['name']

    def __str__(self):
        return str(self.name)


class Project:
    """
    Representing a single project with some of its properties.
    """

    def __init__(self, data: ProjectType, web_host: str):
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
        return [Label(label) for label in self._data['projectLabels']]

    @property
    def count(self) -> Dict[str, int]:
        """
        Project statistics dictionary containing:
        experiments, contributors, children, collaborators, runningExps.
        """
        return self._data['_count']


class Projects:
    """
    Container for a collection of Project objects.
    You can iterate over the projects by for-in loop.
    """

    def __init__(
        self,
        client: Client,
        web_host: str,
        workspace: str,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ):
        self._client = client
        self._web_host = web_host
        self._workspace = workspace
        self._sort = sort
        self._search = search
        self._detail = detail

    def __iter__(self):
        # 按用户遍历情况获取项目信息
        cur_page = 1
        page_size = 20
        while True:
            cur_page += 1
            projects_info = get_entity_projects(
                self._client,
                workspace=self._workspace,
                page=cur_page,
                size=page_size,
                sort=self._sort,
                search=self._search,
                detail=self._detail,
            )
            if cur_page * page_size >= projects_info['total']:
                break

        yield from iter(Project(project, self._web_host) for project in projects_info['list'])
