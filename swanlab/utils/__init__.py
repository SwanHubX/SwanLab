"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13 20:33
@description: 可供外部使用的 SwanLab SDK 工具包
理论上本模块的内容都可以被用户调用，并被写入API文档中
"""

from .experiment import generate_color, generate_id, generate_name
from .metadata import (
    get_metadata_command,
    get_metadata_conda,
    get_metadata_cwd,
    get_metadata_git_branch,
    get_metadata_git_commit,
    get_metadata_git_remote_url,
    get_metadata_hostname,
    get_metadata_os,
    get_metadata_os_pretty,
    get_metadata_pid,
    get_metadata_python_executable,
    get_metadata_python_verbose,
    get_metadata_python_version,
    get_metadata_requirements,
)

__all__ = [
    # experiment
    "generate_color",
    "generate_id",
    "generate_name",
    # metadata
    "get_metadata_os",
    "get_metadata_os_pretty",
    "get_metadata_hostname",
    "get_metadata_pid",
    "get_metadata_cwd",
    "get_metadata_python_version",
    "get_metadata_python_verbose",
    "get_metadata_python_executable",
    "get_metadata_command",
    "get_metadata_git_remote_url",
    "get_metadata_git_branch",
    "get_metadata_git_commit",
    "get_metadata_requirements",
    "get_metadata_conda",
    # monitor
]
