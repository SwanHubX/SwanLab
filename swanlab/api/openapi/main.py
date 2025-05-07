#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""
from swanlab.api.openapi.base import ApiHTTP
from swanlab.api.openapi.experiment import ExperimentAPI
from swanlab.api.openapi.group import GroupAPI
from swanlab.api import code_login
from swanlab.api.http import HTTP
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
            list[dict]: 每个元素是一个字典, 包含工作空间的基础信息:

                - name: str, 工作空间名称
                - username: str, 工作空间唯一标识(用于组织相关的 URL)
                - role: str, 用户在该工作空间中的角色，如 'OWNER' 或 'MEMBER'
        """
        return self.group.list_workspaces()
