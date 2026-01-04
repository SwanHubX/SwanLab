"""
@author: Zhou QiYang
@file: project.py
@time: 2025/12/19 23:49
@description: 定义项目相关的后端API接口
"""

from typing import Optional, List, TYPE_CHECKING

if TYPE_CHECKING:
    from swanlab.core_python.client import Client
from .type import ProjParamType, ProjResponseType


def get_workspace_projects(
    client: "Client",
    *,
    workspace: str,
    page: int = 1,
    size: int = 20,
    sort: Optional[List[str]] = None,
    search: Optional[str] = None,
    detail: Optional[bool] = True,
):
    """
    获取指定页数和条件下的项目信息
    :param client: 已登录的客户端实例
    :param workspace: 工作空间名称
    :param page: 页码
    :param size: 每页项目数量
    :param sort: 排序规则, 可选
    :param search: 搜索的项目名称关键字, 可选
    :param detail: 是否包含项目下实验的相关信息, 可选, 默认为true
    """
    params: ProjParamType = {
        'page': page,
        'size': size,
        'sort': sort,
        'search': search,
        'detail': detail,
    }
    res = client.get(f"/project/{workspace}", params=dict(params))
    projects_info: ProjResponseType = res[0]
    return projects_info
