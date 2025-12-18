"""
@author: Zhou Qiyang
@file: project.py
@time: 2025/12/17 16:24
@description: 项目信息的后端API接口
"""

from typing import Optional, List

from swanlab.core_python.client import Client
from ..type import ProjParamType, ProjResponseType, ProjectType


def get_entity_projects(
    client: Client,
    *,
    username: str,
    sort: Optional[List[str]] = None,
    search: Optional[str] = None,
    detail: Optional[bool] = True,
):
    """
    更新实验状态，注意此接口会将客户端标记为 pending 状态，表示实验已结束
    :param client: 已登录的客户端实例
    :param username: 组织的用户名
    :param sort: 排序规则, 可选
    :param search: 搜索的项目名称关键字, 可选
    :param detail: 是否包含项目下实验的相关信息, 可选, 默认为true
    """
    params: ProjParamType = {
        'page': 1,
        'size': 20,
        'sort': sort,
        'search': search,
        'detail': detail,
    }
    projects: List[ProjectType] = []

    # 循环获取组织中的所有项目
    while True:
        res = client.get(f"/project/{username}", params=dict(params))
        client.pending = True
        projects_info: ProjResponseType = res[0]
        projects.extend(projects_info['list'])
        count = int(projects_info['size']) * int(projects_info['pages'])
        if count >= int(projects_info['total']):
            break
        else:
            params['page'] += 1

    return projects
