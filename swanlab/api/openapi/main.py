#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""

from swanlab.api.openapi.base import ApiHTTP, get_logger
from swanlab.api.openapi.experiment import ExperimentAPI
from swanlab.api.openapi.group import GroupAPI
from swanlab.api import code_login
from swanlab.api.openapi.project import ProjectAPI
from swanlab.api.openapi.types import ApiResponse, Experiment, Pagination, Project
from swanlab.error import KeyFileError
from swanlab.log.log import SwanLog
from swanlab.package import get_key


class OpenApi:
    def __init__(self, key: str = "", log_level: str = "info"):
        self.__logger: SwanLog = get_logger(log_level)

        if key:
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

    def list_workspaces(self) -> ApiResponse[list]:
        """
        获取当前用户的所有工作空间(Group)

        Returns:
            list:
                - code (int): HTTP 错误代码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (list): 一个列表, 其中每个元素是一个字典, 包含相应工作空间的基础信息:
                    - name (str): 工作空间名称
                    - username (str): 工作空间名(用于组织相关的 URL)
                    - role (str): 用户在该工作空间中的角色，如 'OWNER' 或 'MEMBER'
        """
        return self.group.list_workspaces()

    def get_exp_state(self, project: str, exp_cuid: str, username: str = "") -> ApiResponse[dict]:
        """
        获取实验状态

        Args:
           project (str): 项目名
           exp_cuid (str): 实验CUID
           username (str): 工作空间名, 默认为用户个人空间

        Returns:
            dict:
                - code (int): HTTP 错误代码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (dict): 实验状态的字典, 包含以下字段:
                    - state (str): 实验状态, 为 'FINISHED' 或 'RUNNING'
                    - finishedAt (str): 实验完成时间(若有), 格式如 '2024-11-23T12:28:04.286Z'
        """
        return self.experiment.get_exp_state(
            username=username if username else self.http.username, projname=project, expid=exp_cuid  # type: ignore
        )

    def get_experiment(self, project: str, exp_cuid: str, username: str = "") -> ApiResponse[Experiment]:
        """
        获取实验信息

        Args:
            project (str): 项目名
            exp_cuid (str): 实验CUID
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            dict:
                - code (int): HTTP 错误代码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (dict): 实验信息的字典, 包含实验信息
        """
        return self.experiment.get_experiment(
            username=username if username else self.http.username, projname=project, expid=exp_cuid  # type: ignore
        )

    def get_project_exps(
        self, project: str, page: int = 1, size: int = 10, username: str = ""
    ) -> ApiResponse[Pagination[Experiment]]:
        """
        获取项目下的实验列表(分页)

        Args:
            project (str): 项目名
            username (str): 工作空间名, 默认为用户个人空间
            page (int): 页码, 默认为 1
            size (int): 每页大小, 默认为 10

        Returns:
            dict:
                - code (int): HTTP 错误代码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Pagination[Experiment]): 实验列表, 包含以下字段:
                    - total (int): 云端实验总数
                    - list (list[Experiment]): 实验列表, 每个元素都是实验信息的字典, 包含实验信息
                        - 此实验的 profile 只包含 config (实验自定义配置)
        """
        return self.experiment.get_project_exps(
            username=username if username else self.http.username, projname=project, page=page, size=size  # type: ignore
        )

    def list_projects(
        self, username: str = "", detail: bool = True, page: int = 1, size: int = 20
    ) -> ApiResponse[Pagination[Project]]:
        """
        列出一个 workspace 下的所有项目

        Args:
            username (str): 工作空间名, 默认为用户个人空间
            detail (bool): 是否包含项目下实验的相关信息，默认为 True
            page (int): 页码, 默认为 1
            size (int): 每页大小, 默认为 20

        Returns:
            dict: 项目列表信息的字典, 包含以下字段:
                - code (int): HTTP 错误代码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Pagination[Project]): 项目列表, 包含以下字段:
                    - total (int): 指定工作空间下的项目总数
                    - list (list[Project]): 项目列表, 每个元素都是项目信息的字典, 包含项目信息
        """
        username = username or self.http.username  # type: ignore

        return self.project.list_projects(username=username, detail=detail, page=page, size=size)
