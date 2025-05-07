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

    # 获取实验状态，参数对应后台参数命名
    def get_exp_state(self, username: str, projname: str, exp_id: str):
        if not username:
            raise ValueError("The 'username' parameter cannot be empty.")
        if not projname:
            raise ValueError("The 'workspace' parameter cannot be empty.")
        if not exp_id:
            raise ValueError("The 'exp_cuid' parameter cannot be empty.")
        resp = self.http.get(f"/project/{username}/{projname}/runs/{exp_id}/state")
        # 在测试时发现源代码加入了全局异常捕获器，以下代码在评估时请考虑是否添加，给用户做以提示
        if resp is None:
            raise RuntimeError(
                f"Please check username={username},workspace={projname},exp_cuid={exp_id} correctly!"
            )
        return resp
