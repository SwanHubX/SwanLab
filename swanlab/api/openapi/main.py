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
from swanlab.api.openapi.types import ApiResponse, Experiment, Pagination
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
            username=username if username else self.http.username, projname=project, expid=exp_cuid
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
            username=username if username else self.http.username, projname=project, expid=exp_cuid
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
            username=username if username else self.http.username, projname=project, page=page, size=size
        )

    def list_projects(self, username: Optional[str] = None, detail: Optional[bool] = None):
        """
        列出一个 workspace 下的所有项目

        Args:
            username (Optional[str]): 工作空间名, 默认为用户个人空间
            detail (Optional[bool]): 是否包含项目下实验的相关信息，默认为 True

        Returns:
            dict: 项目列表信息的字典, 包含以下字段:

                - total (int): 项目总数
                - list (List[dict]): 项目列表, 每个项目包含以下字段:

                    - cuid (str): 项目的唯一标识符
                    - name (str): 项目名称
                    - description (str): 项目描述
                    - visibility (str): 项目可见性, 如 'PUBLIC' 或 'PRIVATE'
                    - createdAt (str): 项目创建时间, 格式如 '2025-01-06T14:25:29.075Z'
                    - updatedAt (str): 项目更新时间, 格式如 '2025-02-21T09:31:11.473Z'
                    - path (str): 项目路径
                    - group (dict): 项目所属组信息, 包含以下字段:

                        - type (str): 组类型, 如 'PERSON' 或 'TEAM'
                        - username (str): 组用户名
                        - name (str): 组名称(可能为null)

                    - _count (dict): 仅当detail=True时返回, 包含以下字段:

                        - experiments (int): 项目中的实验数量
                        - contributors (int): 项目贡献者数量
                        - children (int): 子项目数量
                        - runningExps (int): 正在运行的实验数量
        """
        username = username or self.http.username
        detail = bool(detail) if detail is not None else True
        if not username:
            return None

        return self.project.list_projects(username=username, detail=detail)
