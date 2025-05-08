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
from swanlab.api.openapi.types import Experiment, ExperimentProfile, ApiErrorResponse


class ExperimentAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, res: dict) -> Experiment:
        return {
            "cuid": res.get("cuid", ""),
            "name": res.get("name", ""),
            "description": res.get("description", ""),
            "state": res.get("state", ""),
            "show": res.get("show", ""),
            "createdAt": res.get("createdAt", ""),
            "finishedAt": res.get("finishedAt", ""),
            "user": {
                "username": res.get("user").get("username", ""),
                "name": res.get("user").get("name", ""),
            },
            "profile": cls.parse_profile(res.get("profile", {})),
        }

    @classmethod
    def parse_profile(cls, res: dict) -> ExperimentProfile:
        return {
            "config": res.get("config", {}),
            "metadata": res.get("metadata", {}),
            "requirements": res.get("requirements", ""),
            "conda": res.get("conda", "")
        }

    def get_exp_state(self, username: str, projname: str, expid: str) -> dict | ApiErrorResponse:
        """
        获取实验状态

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            expid (str): 实验CUID
        """
        return self.http.get(f"/project/{username}/{projname}/runs/{expid}/state")

    def get_experiment(self, username: str, projname: str, expid: str) -> Experiment | ApiErrorResponse:
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
        return ExperimentAPI.parse(resp)

    def get_project_exps(
            self,
            username: str,
            projname: str,
            page: int = 1,
            size: int = 10
    ) -> tuple[list[Experiment], int] | ApiErrorResponse:
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
        return [ExperimentAPI.parse(e) for e in resp.get("list", [])], resp.get("total", 0)
