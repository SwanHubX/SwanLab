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
from swanlab.api.openapi.types import ApiResponse, Experiment, Pagination


class ExperimentAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, body: dict) -> Experiment:
        return Experiment.model_validate({
            "cuid": body.get("cuid") or "",
            "name": body.get("name") or "",
            "description": body.get("description") or "",
            "state": body.get("state") or "",
            "show": body.get("show") or "",
            "createdAt": body.get("createdAt") or "",
            "finishedAt": body.get("finishedAt") or "",
            "user": {
                "username": (body.get("user") or {}).get("username") or "",
                "name": (body.get("user") or {}).get("name") or "",
            },
            "profile": body.get("profile") or {}
        })

    def get_exp_state(self, username: str, projname: str, expid: str) -> ApiResponse[dict]:
        """
        获取实验状态

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            expid (str): 实验CUID
        """
        return self.http.get(f"/project/{username}/{projname}/runs/{expid}/state")

    def get_experiment(self, username: str, projname: str, expid: str) -> ApiResponse[Experiment]:
        """
        获取实验信息

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            expid (str): 实验CUID
        """
        resp = self.http.get(f"/project/{username}/{projname}/runs/{expid}")
        if resp.errmsg:
            return resp
        resp.data = ExperimentAPI.parse(resp.data)
        return resp

    def list_project_exps(
            self,
            username: str,
            projname: str,
            page: int = 1,
            size: int = 10
    ) -> ApiResponse[Pagination[Experiment]]:
        """
        分页获取项目下的实验列表
        该接口返回的实验profile只包含config(用户自定义的配置)

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            page (int): 页码, 默认为1
            size (int): 每页大小, 默认为10
        """
        resp = self.http.get(f"/project/{username}/{projname}/runs", params={"page": page, "size": size})

        if resp.errmsg:
            return resp
        exps = resp.data
        resp.data = Pagination[Experiment].model_validate({
            "total": exps.get("total", 0),
            "list": [ExperimentAPI.parse(e) for e in exps.get("list", [])]
        })
        return resp
