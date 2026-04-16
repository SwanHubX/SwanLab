"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 18:19
@description: SwanLab 运行时 API 封装

绝大多数API使用 Client 对象，少部分API使用requests库直接调用
我们以rpc风格封装API，方便调用
"""

from .experiment import (
    create_or_resume_experiment,
    delete_experiment,
    get_experiment_metrics,
    get_project_experiments,
    get_single_experiment,
    send_experiment_heartbeat,
    update_experiment_state,
)
from .project import delete_project, get_or_create_project, get_project, get_workspace_projects
from .self_hosted import create_user, get_self_hosted_init, get_users
from .user import (
    create_api_key,
    delete_api_key,
    get_api_keys,
    get_latest_api_key,
    get_user_groups,
    get_workspace_info,
)

__all__ = [
    # experiment
    "create_or_resume_experiment",
    "send_experiment_heartbeat",
    "update_experiment_state",
    "get_project_experiments",
    "get_single_experiment",
    "get_experiment_metrics",
    "delete_experiment",
    # project
    "get_project",
    "get_or_create_project",
    "get_workspace_projects",
    "delete_project",
    # user
    "create_api_key",
    "delete_api_key",
    "get_user_groups",
    "get_workspace_info",
    "get_api_keys",
    "get_latest_api_key",
    # self_hosted
    "get_self_hosted_init",
    "create_user",
    "get_users",
]
