"""
@author: Zhou Qiyang
@file: main.py
@time: 2025/12/17 11:39
@description: OpenApi 模块
"""

from typing import Optional, List

from swanlab.api.model import Projects
from swanlab.api.types import ProjectType
from swanlab.core_python import auth, Client
from swanlab.core_python.api.project import get_entity_projects
from swanlab.error import KeyFileError
from swanlab.log import swanlog
from swanlab.package import get_key, HostFormatter

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None


class OpenApi:
    def __init__(self, api_key: Optional[str] = None, host: Optional[str] = None, web_host: Optional[str] = None):
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
            entity: str,
            sort: Optional[List[str]] = None,
            search: Optional[str] = None,
            detail: Optional[bool] = True,
    ) -> Projects:
        projects: List[ProjectType] = get_entity_projects(self._client, username=entity, sort=sort, search=search,
                                                          detail=detail)
        # 扩充一些后端返回中没有的字段
        for proj in projects:
            proj["url"] = f"{self._web_host}/@{proj['path']}"
            print(len(proj))
        result = Projects(projects)

        return result
