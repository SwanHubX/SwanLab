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
        return Experiment.model_validate(
            {
                "cuid": body.get("cuid") or "",
                "name": body.get("name") or "",
                "description": body.get("description") or "",
                "tags": body.get("tags") or "",
                "state": body.get("state") or "",
                "show": bool(body.get("show")),
                "createdAt": body.get("createdAt") or "",
                "finishedAt": body.get("finishedAt") or "",
                "user": {
                    "username": (body.get("user") or {}).get("username") or "",
                    "name": (body.get("user") or {}).get("name") or "",
                },
                "profile": body.get("profile") or {},
            }
        )

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
        resp: ApiResponse = self.http.service.get_exp_info(username=username, projname=projname, expid=expid)
        if resp.errmsg:
            return resp
        resp.data = ExperimentAPI.parse(resp.data)
        return resp

    def list_project_exps(
        self, username: str, projname: str, page: int = 1, size: int = 10
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
        resp.data = Pagination[Experiment].model_validate(
            {"total": exps.get("total", 0), "list": [ExperimentAPI.parse(e) for e in exps.get("list", [])]}
        )
        return resp

    def get_exp_summary(self, expid: str, proid: str, root_expid: str, root_proid: str) -> ApiResponse[dict]:
        """
        获取实验的summary信息
        从House获取, 需要考虑克隆实验

        Args:
            expid (str): 实验CUID
            proid (str): 项目CUID
            root_expid (str): 根实验CUID
            root_proid (str): 根项目CUID
        """
        data = {
            "experimentId": expid,
            "projectId": proid,
        }
        if root_expid and root_proid:
            data["rootExperimentId"] = root_expid
            data["rootProjectId"] = root_proid

        resp = self.http.post(f"/house/metrics/summaries", data=[data])
        if resp.errmsg:
            return resp

        resp.data = list(resp.data.values())[0]
        resp.data = {
            k: {
                "step": v.get("step"),
                "value": v.get("value"),
                "min": {
                    "step": v.get("min").get("index"),
                    "value": v.get("min").get("data"),
                },
                "max": {
                    "step": v.get("max").get("index"),
                    "value": v.get("max").get("data"),
                },
            }
            for k, v in resp.data.items()
        }
        return resp
