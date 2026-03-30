"""
@author: cunyue
@file: git.py
@time: 2026/3/30
@description: Git 信息采集
"""

import subprocess
from typing import Optional


def get_metadata_git_remote_url() -> str:
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


def get_metadata_git_branch() -> str:
    """获取当前分支名"""
    result = subprocess.run(
        ["git", "branch", "--show-current"],
        capture_output=True,
        text=True,
        timeout=5,
        check=True,
    )
    return result.stdout.strip()


def get_metadata_git_commit() -> Optional[str]:
    """获取最新提交 hash"""
    branch = get_metadata_git_branch()
    if not branch:
        return None
    result = subprocess.run(
        ["git", "rev-parse", branch],
        capture_output=True,
        text=True,
        timeout=5,
        check=True,
    )
    return result.stdout.strip()


def parse_git_url(url: str) -> str:
    """将 SSH 格式转换为 HTTPS 格式"""
    if url.startswith("git@"):
        parts = url[4:].split("/", 1)
        host = parts[0].split(":")[0]
        path = parts[1] if len(parts) > 1 else ""
        url = f"https://{host}/{path}"
    return url[:-4] if url.endswith(".git") else url
