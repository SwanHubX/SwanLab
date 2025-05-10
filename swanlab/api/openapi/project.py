#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: project.py
@IDE: pycharm
@Description:
    项目相关的开放API
"""

from swanlab.api.openapi.base import ApiBase, ApiHTTP
from swanlab.api.openapi.types import Project, ApiResponse, Pagination


class ProjectAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, body: dict, detail = True) -> Project:
        project_parser = {
            "cuid": body.get("cuid", ""),
            "name": body.get("name", ""),
            "description": body.get("description", ""),
            "visibility": body.get("visibility", ""),
            "createdAt": body.get("createdAt", ""),
            "updatedAt": body.get("updatedAt", ""),
            "group": {
                "type": body.get("group", {}).get("type", ""),
                "username": body.get("group", {}).get("username", ""),
                "name": body.get("group", {}).get("name", ""),
            },
        }
        if detail:
            project_parser["count"] = body.get("_count")
        return Project.model_validate(project_parser)

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
