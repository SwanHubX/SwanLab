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
        if not projname or not expid:
            raise ValueError("Project name and experiment ID cannot be empty.")
        resp = self.http.get(f"/project/{username}/{projname}/runs/{expid}/state")
        return resp
