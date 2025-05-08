#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: experiment.py
@IDE: pycharm
@Description:
    实验相关的开放API
"""
from swanlab.api.openapi.base import ApiBase, ApiHTTP


class ExperimentAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    def get_exp_state(self, username: str, projname: str, expid: str):
        """
        获取实验状态

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            expid (str): 实验CUID
        """
        return self.http.get(f"/project/{username}/{projname}/runs/{expid}/state")

    def get_experiment(self, username: str, projname: str, expid: str):
        """
        获取实验信息

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            expid (str): 实验CUID
        """
        resp = self.http.get(f"/project/{username}/{projname}/runs/{expid}")
        if "code" in resp:
            return resp
        return {
            "cuid": resp.get("cuid"),
            "name": resp.get("name"),
            "description": resp.get("description"),
            "state": resp.get("state"),
            "show": resp.get("show"),
            "createdAt": resp.get("createdAt"),
            "finishedAt": resp.get("finishedAt"),
            "user": {
                "username": resp.get("user").get("username"),
                "name": resp.get("user").get("name"),
            },
            "profile": resp.get("profile"),
        }

    def get_project_exps(self, username: str, projname: str, page: int = 1, size: int = 10):
        """
        分页获取项目下的实验列表

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            page (int): 页码, 默认为1
            size (int): 每页大小, 默认为10
        """
        resp = self.http.get(f"/project/{username}/{projname}/runs", params={"page": page, "size": size})
        if "code" in resp:
            return resp
        return {
            "total": resp["total"],
            "exps": [
                {
                    "cuid": e.get("cuid"),
                    "name": e.get("name"),
                    "description": e.get("description"),
                    "state": e.get("state"),
                    "show": e.get("show"),
                    "createdAt": e.get("createdAt"),
                    "finishedAt": e.get("finishedAt"),
                    "user": {
                        "username": e.get("user").get("username"),
                        "name": e.get("user").get("name"),
                    },
                    "profile": e.get("profile")
                }
                for e in resp["list"]
            ]
        }