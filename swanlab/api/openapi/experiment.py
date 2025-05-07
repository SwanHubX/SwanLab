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

    # 获取实验状态，参数对应后台参数命名:username为用户名，projname为工作空间的项目名，exp_id为实验id
    def get_exp_state(self, username: str, projname: str, exp_id: str):
        if not projname or not exp_id:
            raise ValueError("'workspace' parameter or 'exp_cuid' parameter is empty")
        else:
            try:
                resp = self.http.get(f"/project/{username}/{projname}/runs/{exp_id}/state")
            except ApiError as e:
                resp = {"code": e.resp.status_code, "message": e.message}
            return resp
