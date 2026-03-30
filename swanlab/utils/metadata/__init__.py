"""
@author: cunyue
@file: __init__.py
@time: 2026/3/30 17:01
@description: SwanLab 元数据相关工具函数和配置
"""

from .conda import get_metadata_conda
from .git import get_metadata_git_branch, get_metadata_git_commit, get_metadata_git_remote_url
from .requirements import get_metadata_requirements
from .runtime import (
    get_metadata_command,
    get_metadata_cwd,
    get_metadata_hostname,
    get_metadata_os,
    get_metadata_os_pretty,
    get_metadata_pid,
    get_metadata_python_executable,
    get_metadata_python_verbose,
    get_metadata_python_version,
)

__all__ = [
    # runtime
    "get_metadata_os",
    "get_metadata_os_pretty",
    "get_metadata_hostname",
    "get_metadata_pid",
    "get_metadata_cwd",
    "get_metadata_python_version",
    "get_metadata_python_verbose",
    "get_metadata_python_executable",
    "get_metadata_command",
    # git
    "get_metadata_git_remote_url",
    "get_metadata_git_branch",
    "get_metadata_git_commit",
    # requirements
    "get_metadata_requirements",
    # conda
    "get_metadata_conda",
]
