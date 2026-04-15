"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 01:18
@description: cmd 模块工具函数
"""

from pathlib import Path

from swanlab.sdk.internal.pkg import fs, nrc
from swanlab.sdk.internal.settings import ROOT_FOLDER, settings
from swanlab.sdk.typings.cmd import LoginType

__all__ = ["get_nrc_path"]


def get_nrc_path(save: LoginType) -> Path:
    """根据登录类型获取 NRC 文件路径"""
    assert save, "LoginType cannot be False when getting NRC path"
    if save == "local":
        return nrc.path(Path.cwd() / ROOT_FOLDER)
    return nrc.path(settings.root)


def append_gitignore(path: Path):
    """在指定路径下创建 .gitignore 文件"""
    if not any(path.iterdir()):
        # 3. 如果为空，写入 .gitignore
        gitignore_path = path / ".gitignore"
        if gitignore_path.exists():
            return
        # 忽略目录下所有文件
        ignore_content = "*"
        fs.safe_write(gitignore_path, ignore_content)
