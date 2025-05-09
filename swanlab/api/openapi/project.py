#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: project.py
@IDE: pycharm
@Description:
    项目相关的开放API
"""

from swanlab.api.http import HTTP
from swanlab.api.openapi.base import ApiBase, ApiHTTP
from swanlab.api.openapi.types import Project, ApiResponse, Pagination


class ProjectAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, body: dict, **kwargs) -> Project:
        project_parser = {
            "cuid": body.get("cuid", ""),
            "name": body.get("name", ""),
            "description": body.get("description", ""),
            "visbility": body.get("visbility", ""),
            "createdAt": body.get("createdAt", ""),
            "updatedAt": body.get("updatedAt", ""),
            "path": body.get("path", ""),
            "group": {
                "type": body.get("group", {}).get("type", ""),
                "username": body.get("group", {}).get("username", ""),
                "name": body.get("group", {}).get("name", ""),
            },
        }
        if kwargs.get("detail", True):
            project_parser["_count"] = body.get("_count")
        return Project.model_validate(project_parser)

    def list_projects(
        self, username: str, detail: bool = True, page: int = 1, size: int = 20
    ) -> ApiResponse[Pagination[Project]]:
        """
        列出一个 workspace 下的所有项目

        Args:
            username: 工作空间名, 默认为用户个人空间
            detail: 是否返回实验详细信息，默认为True
            page: 页码，默认为1
            size: 每页数量，默认为20

        Returns:
            ApiResponse: 包含项目分页数据的响应
        """
        base_url = f"/project/{username}"
        params = {"detail": detail, "page": page, "size": size}

        resp = self.http.get(base_url, params=params)
        if resp.errmsg:
            return resp

        projects = [ProjectAPI.parse(item, detail=detail) for item in resp.data.get("list", [])]
        resp.data = Pagination[Project].model_validate({"total": resp.data.get("total", 0), "list": projects})

        return resp
