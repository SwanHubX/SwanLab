#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/4/29 9:40
@File: __init__.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI包
"""
from typing import Optional, Union, List, Dict

from swanlab.core_python import auth, Client
from swanlab.core_python.api.experiments import get_single_experiment, get_project_experiments
from swanlab.core_python.api.user.self_hosted import get_self_hosted_init
from swanlab.error import KeyFileError, ApiError
from swanlab.log import swanlog
from swanlab.package import HostFormatter, get_key
from .deprecated import OpenApi
from .experiments import Experiments, Experiment
from .projects import Projects
from .user import ApiUser, SuperUser


class Api:
    def __init__(self, api_key: Optional[str] = None, host: Optional[str] = None, web_host: Optional[str] = None):
        """
        初始化 OpenApi 实例，用户需提前登录，或者提供API密钥
        :param api_key: API 密钥，可选
        :param host: API 主机地址，可选
        :param web_host: Web 主机地址，可选
        """
        if host or web_host:
            HostFormatter(host, web_host)()
        if api_key:
            swanlog.debug("Using API key", api_key)
        else:
            swanlog.debug("Using existing key")
            try:
                api_key = get_key()
            except KeyFileError as e:
                swanlog.error("To use SwanLab OpenAPI, please login first.")
                raise RuntimeError("Not logged in.") from e

        self._login_info = auth.code_login(api_key, save_key=False)
        # 一个OpenApi对应一个client，可创建多个api获取从不同的client获取不同账号下的实验信息
        self._client: Client = Client(self._login_info)
        self._web_host = self._login_info.web_host

    def user(self, username: str = None) -> Optional[Union[ApiUser, SuperUser]]:
        # 尝试获取私有化服务信息，如果不是私有化服务，则会报错退出，因为指定user功能仅供私有化用户使用
        try:
            self_hosted_info = get_self_hosted_init(self._client)
        except ApiError as e:
            if username is not None:
                swanlog.error(
                    "You haven't launched a swanlab self-hosted instance. Please check your login status using 'swanlab verify'."
                )
                raise e
            else:
                return ApiUser(self._client, self._login_info)

        if not self_hosted_info["enabled"]:
            raise RuntimeError("SwanLab self-hosted instance hasn't been ready yet.")
        if self_hosted_info["expired"]:
            raise RuntimeError("SwanLab self-hosted instance has expired. Please refresh your licence.")

        # 免费版仅能获取当前api_key登录的用户
        if self_hosted_info["plan"] == 'free':
            if username != self._login_info.username:
                swanlog.warning("Your self-hosted plan is 'free', You will be access to your own account.")
            return ApiUser(self._client, self._login_info)
        # 商业版的根用户可以获取到任何一个用户
        elif self_hosted_info["plan"] == 'commercial':
            if self_hosted_info['root']:
                return SuperUser(self._client, self._login_info, self_hosted=self_hosted_info)
            elif username != self._login_info.username:
                swanlog.warning("Your are not the root user, You will be access to your own account.")
            return ApiUser(self._client, self._login_info)
        # 为教育版预留功能
        else:
            swanlog.warning("The self-hosted plan hasn't been supported yet.")
            return None

    def projects(
        self,
        workspace: str,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ) -> Projects:
        """
        获取指定工作空间（组织）下的所有项目信息
        :param workspace: 工作空间（组织）名称
        :param sort: 排序方式，可选
        :param search: 搜索关键词，可选
        :param detail: 是否返回详细信息，可选
        :return: Projects 实例，可遍历获取项目信息
        """
        return Projects(
            client=self._client,
            web_host=self._web_host,
            workspace=workspace,
            sort=sort,
            search=search,
            detail=detail,
        )

    def runs(self, path: str, filters: Dict[str, object] = None) -> Experiments:
        """
        获取指定项目下的所有实验信息
        :param path: 项目路径，格式为 'username/project'
        :return: Experiments 实例，可遍历获取实验信息
        :param filters: 筛选实验的条件，可选
        """
        return Experiments(client=self._client, path=path, web_host=self._web_host, filters=filters)

    def run(
        self,
        path: str,
    ) -> Experiment:
        """
        获取指定实验的信息
        :param path: 实验路径，格式为 'username/project/run_id'
        :return: Experiment 实例，包含实验信息
        """
        # TODO: 待后端完善后替换成专用的接口
        if len(path.split('/')) != 3:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project/run_id'")
        _data = get_single_experiment(self._client, path=path)
        proj_path = path.rsplit('/', 1)[0]
        data = get_project_experiments(
            self._client, path=proj_path, filters={'name': _data['name'], 'created_at': _data['createdAt']}
        )
        return Experiment(data=data[0], client=self._client, path=proj_path, web_host=self._web_host, line_count=1)


__all__ = ["Api", "OpenApi"]
