"""
@author: Zhou QiYang
@file: __init__.py
@time: 2026/1/11 16:31
@description: OpenApi 中的项目对象迭代器
"""

from typing import List, Optional, Iterator

from swanlab.api.project import Project
from swanlab.core_python.api.project import get_workspace_projects
from swanlab.core_python.api.type import ProjResponseType
from swanlab.core_python.client import Client


class Projects:
    """
    Container for a collection of Project objects.
    You can iterate over the projects by for-in loop.
    """

    def __init__(
        self,
        client: Client,
        *,
        web_host: str,
        workspace: str,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ) -> None:
        self._client = client
        self._web_host = web_host
        self._workspace = workspace
        self._sort = sort
        self._search = search
        self._detail = detail

    def __iter__(self) -> Iterator[Project]:
        # 按用户遍历情况获取项目信息
        cur_page = 0
        while True:
            cur_page += 1
            resp: ProjResponseType = get_workspace_projects(
                self._client,
                workspace=self._workspace,
                page=cur_page,
                size=20,
                sort=self._sort,
                search=self._search,
                detail=self._detail,
            )
            for p in resp['list']:
                yield Project(self._client, data=p, web_host=self._web_host)

            if cur_page >= resp['pages']:
                break


__all__ = ["Projects"]
