#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""
from typing import Dict, List

from swanlab.api import code_login
from swanlab.api.openapi.base import ApiHTTP, fetch_paginated_api, get_logger
from swanlab.api.openapi.experiment import ExperimentAPI
from swanlab.api.openapi.group import GroupAPI
from swanlab.api.openapi.project import ProjectAPI
from swanlab.api.openapi.types import ApiResponse, Experiment, Project
from swanlab.error import KeyFileError
from swanlab.log.log import SwanLog
from swanlab.package import get_key


class OpenApi:
    def __init__(self, api_key: str = "", log_level: str = "info"):
        self.__logger: SwanLog = get_logger(log_level)

        if api_key:
            self.__logger.debug("Using API key", api_key)
            self.__key = api_key
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

    def list_workspaces(self) -> ApiResponse[List[Dict]]:
        """
        获取当前用户的所有工作空间(Group)

        Returns:
            ApiResponse[List]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (List[Dict]): 一个列表, 其中每个元素是一个字典, 包含相应工作空间的基础信息:
                    - name (str): 工作空间名称
                    - username (str): 工作空间名(用于组织相关的 URL)
                    - role (str): 用户在该工作空间中的角色，如 'OWNER' 或 'MEMBER'
        """
        return self.group.list_workspaces()

    def get_exp_state(
            self, project: str, exp_cuid: str, username: str = ""
    ) -> ApiResponse[Dict]:
        """
        获取实验状态

        Args:
           project (str): 项目名
           exp_cuid (str): 实验CUID
           username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[Dict]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Dict): 实验状态的字典, 包含以下字段:
                    - state (str): 实验状态, 为 'FINISHED' 或 'RUNNING'
                    - finishedAt (str): 实验完成时间(若有), 格式如 '2024-11-23T12:28:04.286Z'
        """
        return self.experiment.get_exp_state(
            username=username if username else self.http.username, projname=project, expid=exp_cuid
        )

    def get_experiment(
            self, project: str, exp_cuid: str, username: str = ""
    ) -> ApiResponse[Experiment]:
        """
        获取实验信息

        Args:
            project (str): 项目名
            exp_cuid (str): 实验CUID
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[Experiment]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Dict): 实验信息的字典, 包含实验信息
        """
        return self.experiment.get_experiment(
            username=username if username else self.http.username, projname=project, expid=exp_cuid
        )

    def list_project_exps(
        self, project: str, username: str = ""
    ) -> ApiResponse[List[Experiment]]:
        """
        获取一个项目下的所有实验

        Args:
            project (str): 项目名
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[Experiment]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (List[Experiment]): 实验列表, 每个元素包含一个实验的信息
                    - 此实验的 profile 只包含 config (实验自定义配置)
        """
        return fetch_paginated_api(
            api_func=self.experiment.list_project_exps,
            projname=project,
            username=username if username else self.http.username
        )

    def list_projects(
        self, username: str = "", detail: bool = True
    ) -> ApiResponse[List[Project]]:
        """
        获取一个工作空间下的所有项目

        Args:
            username (str): 工作空间名, 默认为用户个人空间
            detail (bool): 是否包含详细统计信息，默认为 True

        Returns:
            ApiResponse[List[Project]]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (List[Project]): 项目列表, 每个元素包含一个项目的信息
        """
        return fetch_paginated_api(
            api_func=self.project.list_projects,
            username=username if username else self.http.username,
            detail=detail
        )
