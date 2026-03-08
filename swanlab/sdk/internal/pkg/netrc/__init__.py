"""
@author: cunyue
@file: __init__.py
@time: 2026/3/8 15:42
@description: SwanLab SDK Netrc 工具，提供纯粹的凭证底层文件读写
"""

import netrc
import os
import stat
from pathlib import Path
from typing import Optional, Tuple

__all__ = ["get_nrc_path", "remove_host_suffix", "read_netrc_by_host", "write_netrc"]


def get_nrc_path(path: Path) -> Path:
    """获取netrc文件路径，并不保证文件、目录存在"""
    return path / ".netrc"


def remove_host_suffix(host: str, *suffixes: str) -> str:
    """移除host的后缀"""
    host = host.rstrip()
    if not suffixes:
        return host
    for suffix in suffixes:
        if not suffix:
            continue
        if host.endswith(suffix):
            return host[: -len(suffix)]
    return host


def read_netrc_by_host(nrc_path: Path, target_host: str) -> Optional[Tuple[str, str]]:
    """
    读取 netrc 中指定 host 的身份信息。
    :return: 成功时返回 (username, password/api_key)，失败返回 None
    """
    if not nrc_path.exists():
        return None

    try:
        nrc = netrc.netrc(nrc_path)
    # 增加 PermissionError 以兼容 Windows 这种“把目录当文件开”会报权限错误的特性
    except (netrc.NetrcParseError, IsADirectoryError, PermissionError):
        raise IOError(
            f"Failed to access or parse netrc file at {nrc_path}. "
            f"Please check if the path is a directory, has incorrect permissions, or contains syntax errors."
        )

    info = nrc.authenticators(target_host)

    # 向下兼容逻辑：处理旧版本写入了带 /api 后缀的情况
    if info is None:
        legacy_host = target_host + "/api"
        info = nrc.authenticators(legacy_host)
        if info is not None:
            # 根据全局单点登录原则，直接全量覆盖刷新文件
            nrc.hosts = {target_host: info}

            with open(nrc_path, "w", encoding="utf-8") as f:
                f.write(repr(nrc))

    if info is None:
        return None

    # info 的结构是 (login, account, password)
    return info[0], info[2]


def write_netrc(nrc_path: Path, host: str, username: str, password: str) -> None:
    """
    将凭证写入 netrc 文件。
    遵循全局单点登录原则：每次写入都会清空其他所有历史/不同环境的登录凭证。
    """
    # 1. 核心防御：处理路径冲突
    if nrc_path.exists() and nrc_path.is_dir():
        raise IOError(
            f"The path {nrc_path} is a directory, but a file is expected. "
            f"Please remove the directory or change the SWANLAB_ROOT settings."
        )

    # 2. 确保父目录存在
    nrc_path.parent.mkdir(parents=True, exist_ok=True)

    # 3. 文件初始化与权限控制 (仅属主读写)
    if not nrc_path.exists():
        try:
            nrc_path.touch(mode=0o600)
        except OSError as e:
            raise OSError(f"Failed to create netrc file at {nrc_path}: {e}. Check your permissions.")
    else:
        os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)

    # 4. 尝试加载或重置 netrc 对象
    try:
        nrc = netrc.netrc(nrc_path)
    except (netrc.NetrcParseError, IOError):
        nrc_path.write_text("")
        nrc = netrc.netrc(nrc_path)

    # 5. 写入逻辑：全局单点登录，直接覆写整个 hosts 字典
    new_info = (username, "", password)
    info = nrc.authenticators(host)

    if info is None or (info[0], info[2]) != (new_info[0], new_info[2]):
        nrc.hosts = {host: new_info}  # <-- 强制覆写，顶掉旧账号
        try:
            with open(nrc_path, "w", encoding="utf-8") as f:
                f.write(repr(nrc))
            os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)
        except OSError as e:
            raise OSError(f"Failed to write API Key to {nrc_path}: {e}")
