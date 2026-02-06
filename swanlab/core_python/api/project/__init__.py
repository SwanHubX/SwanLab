"""
@author: Zhou QiYang
@file: __init__.py
@time: 2025/12/19 23:49
@description: 定义项目相关的后端API接口
"""

from typing import Optional, List, TYPE_CHECKING

from swanlab.core_python.api.type import ProjResponseType, ProjectType

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


def get_workspace_projects(
    client: "Client",
    *,
    path: str,
    page: int = 1,
    size: int = 20,
    sort: Optional[str] = None,
    search: Optional[str] = None,
    detail: Optional[bool] = True,
) -> ProjResponseType:
    """
    获取指定页数和条件下的项目信息
    :param client: 已登录的客户端实例
    :param path: 工作空间名称
    :param page: 页码
    :param size: 每页项目数量
    :param sort: 排序规则, 可选
    :param search: 搜索的项目名称关键字, 可选
    :param detail: 是否包含项目下实验的相关信息, 可选, 默认为true
    """
    params = {
        'page': page,
        'size': size,
        'sort': sort,
        'search': search,
        'detail': detail,
    }
    res = client.get(f"/project/{path}", params=dict(params))
    return res[0]


def get_project_info(client: "Client", *, path: str) -> ProjectType:
    """
    获取指定路径的项目信息
    :param client: 已登录的客户端实例
    :param path: 项目路径 'username/project'
    """
    res = client.get(f"/project/{path}")
    return res[0]


__all__ = ["get_workspace_projects", "get_project_info"]
