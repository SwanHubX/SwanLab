#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/4/30 12:11
@File: group.py
@IDE: pycharm
@Description:
    组织相关的开放API
"""

from swanlab.api.openapi.base import ApiBase, ApiHTTP
from swanlab.api.openapi.types import ApiResponse


class GroupAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    def list_workspaces(self) -> ApiResponse[list]:
        resp = self.http.get("/group/")
        if resp.errmsg:
            return resp
        groups = resp.data.get("list", [])
        resp.data = [
            {
                "name": item["name"],
                "username": item["username"],
                "role": item["role"]
            }
            for item in groups
        ]
        return resp
