"""
@author: cunyue
@file: git.py
@time: 2026/3/30
@description: Git 信息采集
"""

import subprocess
from typing import Optional

from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.internal.probe_python.typings import GitSnapshot


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
    """将 Git 远程地址归一化为展示用 HTTPS 形式，并去掉 ``.git`` 后缀。

    - SSH scp-like: ``git@host:owner/repo[.git]`` → ``https://host/owner/repo``
    - SSH 带端口:   ``git@host:port/owner/repo[.git]`` → ``https://host:port/owner/repo``

      host 段冒号右侧为纯数字时视为端口予以保留，否则当作路径首段（标准 scp-like 写法）。
    - 其他 (HTTPS / ssh:// / git:// ...) 原样返回。
    - 统一去掉末尾 ``.git`` 后缀。
    """
    if url.startswith("git@"):
        # 去掉 'git@' 后按首个 '/' 切分: 之前是 host 段, 之后是仓库路径
        parts = url[4:].split("/", 1)
        host, path = parts[0], parts[1] if len(parts) > 1 else ""
        if ":" in host:
            # host 段含冒号: 右侧纯数字 → 视为 host:port 保留; 否则当作路径首段
            host, tail = host.rsplit(":", 1)
            if tail.isdigit():
                host = f"{host}:{tail}"
            else:
                path = f"{tail}/{path}" if path else tail
        # path 为空时不拼尾部斜杠, 避免生成 "https://host/" 这类畸形结果
        url = f"https://{host}/{path}" if path else f"https://{host}"
    if url.endswith(".git"):
        url = url[:-4]
    return url
