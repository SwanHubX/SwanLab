#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""
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
        self.http: HTTP = HTTP(self.login_info)

        self.group = GroupAPI(self.http)
        self.experiment = ExperimentAPI(self.http)
        self.project = ProjectAPI(self.http)

    def list_workspaces(self):
        return self.group.list_workspaces()
