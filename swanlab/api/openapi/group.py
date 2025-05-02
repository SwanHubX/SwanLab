#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/4/30 12:11
@File: group.py
@IDE: pycharm
@Description:
    组织相关的开放API
"""

from swanlab.api.http import HTTP
from swanlab.api.openapi.base import ApiBase


class GroupAPI(ApiBase):
    def __init__(self, http: HTTP):
        super().__init__(http)

    def list_workspaces(self):
        """
        获取当前用户的所有工作空间(Group)

        :return: 一个列表，每个元素是一个字典，包含工作空间的基础信息：
            [
                {
                    "name": str,       # 工作空间名称
                    "username": str,   # 工作空间唯一标识(用于组织相关的URL)
                    "role": str        # 用户在该工作空间中的角色，例如 'OWNER' 或 'MEMBER'
                },
                ...
            ]
        """
        resp = self.http.get("/group/")
        groups: list = [
            {
                "name": item["name"],
                "username": item["username"],
                "role": item["role"]
            }
            for item in resp.get("list", [])
        ]
        return groups
