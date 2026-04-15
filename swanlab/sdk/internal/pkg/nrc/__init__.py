"""
@author: cunyue
@file: __init__.py
@time: 2026/4/15 21:16
@description: Netrc 工具，提供纯粹的凭证底层文件读写
"""

import netrc
import os
import stat
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import urlsplit

from swanlab.sdk.internal.pkg import safe

__all__ = ["fmt", "read", "write"]


def path(p: Path):
    """
    获取 netrc 文件的完整路径
    """
    return p / ".netrc"


def fmt(host: str) -> str:
    """
    格式化域名，移除路由等后缀，仅保留协议、端口、域名，如果没有协议，自动添加 https 协议

    :param host: 需要格式化的 host 字符串

    :raises ValueError: 当 host 为空或仅包含空白字符时抛出
    """
    host = host.strip()
    if not host:
        raise ValueError("Host cannot be empty or whitespace.")

    has_scheme = "://" in host
    parsed = urlsplit(host if has_scheme else f"//{host}")

    if not parsed.hostname:
        return host

    scheme = parsed.scheme if has_scheme and parsed.scheme else "https"
    result = f"{scheme}://{parsed.hostname}"

    if parsed.port is not None:
        result += f":{parsed.port}"

    return result


def read(path: Path) -> Optional[Tuple[str, str, str]]:
    """
    读取 netrc 文件中的凭证信息，依次返回 api_key, api_host, web_host

    映射规则（与写入一致）：
        - machine (host) -> api_host
        - login (username) -> web_host
        - password -> api_key

    :param path: netrc 文件路径
    :return: (api_key, api_host, web_host)，读取失败返回 None
    """

    if not path.exists():
        return None

    try:
        nrc = netrc.netrc(path)
    except (netrc.NetrcParseError, IsADirectoryError, PermissionError):
        raise IOError(
            f"Failed to access or parse netrc file at {path}. "
            f"Please check if the path is a directory, has incorrect permissions, or contains syntax errors."
        )

    if not nrc.hosts:
        return None

    with safe.block(message="Failed to read credentials from netrc file"):
        # 由于"全局单点登录"设计，.netrc 中理论上只有一个 host，直接取第一个
        machine = list(nrc.hosts.keys())[0]
        login, _, password = nrc.hosts[machine]

        # 向下兼容：处理旧版本写入了带 /api 后缀的 machine
        if machine.endswith("/api"):
            machine = machine[: -len("/api")]
            # 自愈：将修正后的 host 写回文件
            nrc.hosts = {machine: (login, _, password)}
            with open(path, "w", encoding="utf-8") as f:
                f.write(repr(nrc))

        api_host = fmt(machine)
        web_host = fmt(login) if login else ""
        return password, api_host, web_host

    return None


def write(nrc_path: Path, host: str, username: str, password: str):
    """
    将凭证写入 netrc 文件。
    遵循全局单点登录原则：每次写入都会清空其他所有历史/不同环境的登录凭证。
    """
    # 1. 核心防御
    # 检查字段是否合法
    for field_name, value in [("host", host), ("username", username), ("password", password)]:
        if any(c in value for c in ("\n", "\r", " ", "\t")):
            raise ValueError(f"Invalid characters in {field_name}. Newlines and spaces are not allowed.")
    # 核心防御：处理路径冲突
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
        with safe.block(message="Failed to write credentials to netrc file"):
            with open(nrc_path, "w", encoding="utf-8") as f:
                f.write(repr(nrc))
            os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)


def remove(nrc_path: Path) -> None:
    """
    删除 netrc 文件中的所有凭证条目。
    遵循全局单点登录原则：直接清空整个 hosts 字典并写回文件。
    如果文件不存在则不做任何操作。
    """
    if not nrc_path.exists():
        return

    try:
        nrc_obj = netrc.netrc(nrc_path)
    except (netrc.NetrcParseError, IsADirectoryError, PermissionError):
        # 文件损坏或不可读，直接删除文件
        nrc_path.unlink(missing_ok=True)
        return

    if not nrc_obj.hosts:
        return

    nrc_obj.hosts = {}
    with open(nrc_path, "w", encoding="utf-8") as f:
        f.write(repr(nrc_obj))
    os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)
