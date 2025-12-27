"""
@author: Zhou Qiyang
@file: main.py
@time: 2025/12/17 11:39
@description: OpenApi 模块
"""

from typing import Optional, List

from swanlab.core_python import auth, Client
from swanlab.error import KeyFileError
from swanlab.log import swanlog
from swanlab.package import get_key, HostFormatter
from swanlab.core_python.api.experiment import get_single_experiment
from .model import Projects, Experiments, Experiment

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None


class OpenApi:
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

        login_info = auth.code_login(api_key, save_key=False)
        # 一个OpenApi对应一个client，可创建多个api获取从不同的client获取不同账号下的实验信息
        self._client: Client = Client(login_info)
        self._web_host = login_info.web_host

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

    def runs(
        self,
        path: str,
    ) -> Experiments:
        """
        获取指定项目下的所有实验信息
        :param path: 项目路径，格式为 'username/project'
        :return: Experiments 实例，可遍历获取实验信息
        """
        return Experiments(
            client=self._client,
            path=path,
            web_host=self._web_host,
        )

    def run(
        self,
        path: str,
    ) -> Experiment:
        """
        获取指定实验的信息
        :param path: 实验路径，格式为 'username/project/expid'
        :return: Experiment 实例，包含实验信息
        """
        if len(path.split('/')) != 3:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project/expid'")
        data = get_single_experiment(self._client, path=path)
        return Experiment(data=data, path=path.rsplit('/', 1)[0], web_host=self._web_host, line_count=1)
