#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: project.py
@IDE: pycharm
@Description:
    项目相关的开放API
"""

from swanlab.api.base import ApiBase, ApiHTTP
from swanlab.api.types import ApiResponse, Pagination, Project


class ProjectAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, body: dict, detail=True) -> Project:
        project_parser = {
            "cuid": body.get("cuid") or "",
            "name": body.get("name") or "",
            "description": body.get("description") or "",
            "visibility": body.get("visibility") or "",
            "createdAt": body.get("createdAt") or "",
            "updatedAt": body.get("updatedAt") or "",
            "group": {
                "type": (body.get("group") or {}).get("type") or "",
                "username": (body.get("group") or {}).get("username") or "",
                "name": (body.get("group") or {}).get("name") or "",
            },
        }
        if detail:
            project_parser["count"] = body.get("_count") or {}
        return Project.model_validate(project_parser)

    def delete_project(
        self, username: str, project: str
    ) -> ApiResponse[None]:
        """
        删除一个项目

        Args:
            username (str): 工作空间名
            project (str): 项目名
        """
        return self.http.delete(f"/project/{username}/{project}", params={})

    def list_projects(
        self, username: str, detail = True, page: int = 1, size: int = 10
    ) -> ApiResponse[Pagination[Project]]:
        """
        列出一个 workspace 下的所有项目

        Args:
            username (str): 工作空间名, 默认为用户个人空间
            detail (bool): 是否返回详细统计信息，默认为 True
            page (int): 页码，默认为 1
            size (int): 每页数量，默认为 10
        """
        resp = self.http.get(f"/project/{username}", params={"detail": detail, "page": page, "size": size})
        if resp.errmsg:
            return resp

        resp.data = Pagination[Project].model_validate({
            "total": resp.data.get("total", 0),
            "list": [ProjectAPI.parse(project, detail=detail) for project in resp.data.get("list", [])]
        })

        return resp
