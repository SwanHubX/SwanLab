#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""
from typing import Optional

from swanlab.api.openapi.base import ApiHTTP
from swanlab.api.openapi.experiment import ExperimentAPI
from swanlab.api.openapi.group import GroupAPI
from swanlab.api import code_login
from swanlab.api.openapi.project import ProjectAPI
from swanlab.error import KeyFileError
from swanlab.log.log import SwanLog
from swanlab.package import get_key


class OpenApi:
    def __init__(self, key: str = None, log_level: str = "info"):
        self.__logger: SwanLog = SwanLog("swanlab.openapi", log_level)

        if key is not None:
            self.__logger.debug("Using API key", key)
            self.__key = key
            self.login_info = code_login(self.__key, False)
        else:
            self.__logger.debug("Using existing key")
            try:
                self.__key = get_key()
            except KeyFileError as e:
                self.__logger.error("To use SwanLab OpenAPI, please login first.")
                raise RuntimeError("Not logged in.") from e
            self.login_info = code_login(self.__key, False)

        self.username = self.login_info.username
        self.__http: ApiHTTP = ApiHTTP(self.login_info)

        self.group = GroupAPI(self.http)
        self.experiment = ExperimentAPI(self.http)
        self.project = ProjectAPI(self.http)

    @property
    def http(self) -> ApiHTTP:
        """
        当前使用的ApiHTTP对象
        """
        return self.__http

    def list_workspaces(self):
        """
        获取当前用户的所有工作空间(Group)

        Returns:
            list[dict]: 一个列表, 其中每个元素是一个字典, 包含相应工作空间的基础信息:

                - name (str): 工作空间名称
                - username (str): 工作空间名(用于组织相关的 URL)
                - role (str): 用户在该工作空间中的角色，如 'OWNER' 或 'MEMBER'

            若请求失败, 将返回包含以下字段的字典:

                - code (int): HTTP 错误代码
                - message (str): 错误信息
        """
        return self.group.list_workspaces()

    def get_exp_state(self, project: str, exp_cuid: str, username: Optional[str] = None):
        """
        获取实验状态

        Args:
           project (str): 项目名
           exp_cuid (str): 实验id
           username (Optional[str]): 工作空间名, 默认为用户个人空间

        Returns:
            dict: 实验状态的字典, 包含以下字段:

                - state (str): 实验状态, 为 'FINISHED' 或 'RUNNING'
                - finishedAt (str): 实验完成时间（若有）, 格式如 '2024-11-23T12:28:04.286Z'

            若请求失败, 将返回包含以下字段的字典:

                - code (int): HTTP 错误代码
                - message (str): 错误信息
        """
        if username:
            return self.experiment.get_exp_state(username=username, projname=project, expid=exp_cuid)
        else:
            return self.experiment.get_exp_state(username=self.http.username, projname=project, expid=exp_cuid)
