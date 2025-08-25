#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:00
@File: main.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""
from typing import Dict, List, Union

from swanlab.api.base import ApiHTTP, get_logger
from swanlab.api.experiment import ExperimentAPI
from swanlab.api.group import GroupAPI
from swanlab.api.project import ProjectAPI
from swanlab.api.types import ApiResponse, Experiment, Project
from swanlab.core_python import auth
from swanlab.error import KeyFileError
from swanlab.log.log import SwanLog
from swanlab.package import get_key

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None


class OpenApi:
    def __init__(self, api_key: str = "", log_level: str = "info"):
        self.__logger: SwanLog = get_logger(log_level)

        if api_key:
            self.__logger.debug("Using API key", api_key)
            self.__key = api_key
            self.login_info = auth.code_login(self.__key, False)
        else:
            self.__logger.debug("Using existing key")
            try:
                self.__key = get_key()
            except KeyFileError as e:
                self.__logger.error("To use SwanLab OpenAPI, please login first.")
                raise RuntimeError("Not logged in.") from e
            self.login_info = auth.code_login(self.__key, False)

        self.username = self.login_info.username
        self.http: ApiHTTP = ApiHTTP(self.login_info)
        self.service = self.http.service

        self.group = GroupAPI(self.http)
        self.experiment = ExperimentAPI(self.http)
        self.project = ProjectAPI(self.http)

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

    def get_experiment(
            self,
            project: str,
            exp_id: str,
            username: str = ""
    ) -> ApiResponse[Experiment]:
        """
        获取实验信息

        Args:
            project (str): 项目名
            exp_id (str): 实验CUID
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[Experiment]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Dict): 实验信息的字典, 包含实验信息
        """
        return self.experiment.get_experiment(
            username=username if username else self.http.username, projname=project, exp_id=exp_id
        )

    def delete_experiment(
            self,
            project: str,
            exp_id: str,
            username: str = ""
    ) -> ApiResponse[None]:
        """
        删除实验

        Args:
            project (str): 项目名
            exp_id (str): 实验CUID
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[None]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (None): 无数据返回
        """
        return self.experiment.delete_experiment(
            username=username if username else self.http.username, projname=project, exp_id=exp_id
        )

    def list_experiments(
            self,
            project: str,
            username: str = ""
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
        return self.service.fetch_paginated_api(
            api_func=self.experiment.list_experiments,
            projname=project,
            username=username if username else self.http.username
        )

    def delete_project(
            self,
            project: str,
            username: str = ""
    ) -> ApiResponse[None]:
        """
        删除一个项目

        Args:
            project (str): 项目名
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[None]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (None): 无数据返回
        """
        return self.project.delete_project(
            username=username if username else self.http.username, project=project
        )

    def list_projects(
            self,
            username: str = "",
            detail: bool = True
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
        return self.service.fetch_paginated_api(
            api_func=self.project.list_projects,
            username=username if username else self.http.username,
            detail=detail
        )

    def get_summary(
            self,
            project: str,
            exp_id: str,
            username: str = ""
    ) -> ApiResponse[Dict]:
        """
        获取实验的概要信息

        Args:
            project (str): 项目名
            exp_id (str): 实验CUID
            username (str): 工作空间名, 默认为用户个人空间

        Returns:
            ApiResponse[Dict]:
                - code (int): HTTP 状态码
                - errmsg (str): 错误信息, 仅在请求有错误时非空
                - data (Dict): 实验的概要信息字典, 包含用户训练各指标的最大最小值, 及其对应步数
        """
        username = username if username else self.http.username
        project_cuid = self.service.get_project_info(username=username, projname=project).data.get("cuid", "")
        exp = self.service.get_exp_info(username=username, project=project, exp_id=exp_id)
        return self.experiment.get_summary(
            exp_id=exp_id,
            pro_id=project_cuid,
            root_exp_id=exp.data.get("rootExpId", ""), 
            root_pro_id=exp.data.get("rootProId", "")
        )

    def get_metrics(
        self,
        exp_id: str,
        keys: Union[str, List[str]],
    ) -> ApiResponse[DataFrame]:
        """
        获取实验的指标数据

        Args:
            exp_id (str): 实验CUID
            keys (str | List[str]): 指标key, 单个字符串或字符串列表

        Returns:
            ApiResponse[DataFrame]: 包含指标数据的响应, 指标数据以 DataFrame 格式返回
            在DataFrame中, 每个key对应两个列, 分别为key和key_timestamp, 表示指标值和时间戳
        """
        if isinstance(keys, str):
            keys = [keys]
        return self.experiment.get_metrics(exp_id, keys)

    def get_exp_summary(self, *args, **kwargs) -> ApiResponse[Dict]:
        """
        获取实验的概要信息
        @deprecated, 请使用 `get_experiment_summary`

        Returns:
            ApiResponse[Dict]: 包含实验概要信息的响应
        """
        return self.get_summary(*args, **kwargs)

    def list_project_exps(self, *args, **kwargs) -> ApiResponse[List[Experiment]]:
        """
        获取一个项目下的所有实验
        @deprecated, 请使用 `list_experiments`

        Returns:
            ApiResponse[List[Experiment]]: 包含实验列表的响应
        """
        return self.list_experiments(*args, **kwargs)
