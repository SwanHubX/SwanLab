"""
@author: cunyue
@file: git.py
@time: 2026/3/30
@description: Git 信息采集
"""

import subprocess
from pathlib import PurePosixPath
from typing import Optional

from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.run.system import GitSnapshot


def get() -> GitSnapshot:
    """创建 GitSnapshot 对象"""
    return GitSnapshot(
        remote_url=get_remote_url(),
        branch=get_branch(),
        commit=get_commit(),
    )


@safe.decorator(level="debug", message="Failed to get git remote url")
def get_remote_url() -> str:
    """获取 Git 远程仓库地址"""
    result = subprocess.run(
        ["git", "config", "--get", "remote.origin.url"],
        capture_output=True,
        text=True,
        timeout=5,
        check=True,
    )
    url = result.stdout.strip()
    return parse_git_url(url)


@safe.decorator(level="debug", message="Failed to get git branch")
def get_branch() -> str:
    """获取当前分支名"""
    result = subprocess.run(
        ["git", "branch", "--show-current"],
        capture_output=True,
        text=True,
        timeout=5,
        check=True,
    )
    return result.stdout.strip()


@safe.decorator(level="debug", message="Failed to get git commit")
def get_commit() -> Optional[str]:
    """获取最新提交 hash"""
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        timeout=5,
        check=True,
    )
    return result.stdout.strip()


def parse_git_url(url: str) -> str:
    """将 SSH 格式转换为 HTTPS 格式"""
    if not url.startswith("git@"):
        return url

    # 1. 去掉 'git@' 协议头
    # 2. 将第一个 ':' 替换为 '/'，使其符合路径规范
    # 比如: github.com:SwanHubX/SwanLab.git -> github.com/SwanHubX/SwanLab.git
    normalized_path = url[4:].replace(":", "/", 1)

    # 3. 使用 PurePosixPath 包装 (强制使用正斜杠)
    path_obj = PurePosixPath(normalized_path)

    # 4. 重新组装成 HTTPS URL
    # path_obj 此时代表了 github.com/SwanHubX/SwanLab.git
    return f"https://{path_obj}"
