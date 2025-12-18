"""
@author: Zhou Qiyang
@file: model.py
@time: 2025/12/18 20:10
@description: OpenApi查询结果将以对象返回，并且对后端的返回字段进行一些筛选
"""

from typing import List, Dict

from .type import ProjectType, ProjectLabelType, GroupType


class Group:
    def __init__(self, data: GroupType):
        self._data = data

    @property
    def username(self):
        return self._data['username']

    @property
    def status(self):
        return self._data['status']

    @property
    def type(self):
        return self._data['type']


class ProjectLabel:
    def __init__(self, data: ProjectLabelType):
        self._data = data

    @property
    def name(self):
        return self._data['name']


class Project:
    def __init__(self, data: ProjectType):
        self._data = data

    @property
    def cuid(self):
        return self._data['cuid']

    @property
    def name(self):
        return self._data['name']

    @property
    def path(self):
        return self._data['path']

    @property
    def url(self):
        return self._data['url']

    @property
    def description(self):
        return self._data['description']

    @property
    def visibility(self):
        return self._data['visibility']

    @property
    def createdAt(self):
        return self._data['createdAt']

    @property
    def updatedAt(self):
        return self._data['updatedAt']

    @property
    def projectLabels(self) -> List[ProjectLabel]:
        return [ProjectLabel(label) for label in self._data['projectLabels']]

    @property
    def group(self) -> Group:
        return Group(self._data['group'])

    @property
    def count(self) -> Dict[str, int]:
        return self._data['_count']


class Projects:
    def __init__(self, data: List[ProjectType]):
        self._projects: List[Project] = [Project(d) for d in data]

    def __iter__(self):
        return iter(self._projects)
