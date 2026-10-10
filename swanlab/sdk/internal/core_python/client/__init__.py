"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:04
@description: SwanLab 运行时客户端，用于与 SwanLab 后端 API 进行交互
"""

import threading
from typing import Optional

from swanlab.sdk.internal.pkg import client, console, fork
from swanlab.sdk.typings.pkg.client import JSONBody, JSONDict

__all__ = ["exists", "reset", "new", "get", "post", "put", "patch", "delete"]

# ==============================================================================
# 模块级全局状态与代理快捷函数
# ==============================================================================

_default_client: Optional[client.Client] = None
# 共享锁：保护 client 单例的创建与销毁，覆盖网络认证构造期间
_client_lock = threading.Lock()


def _after_fork() -> None:
    global _client_lock, _default_client
    _client_lock = threading.Lock()
    _default_client = None


fork.register(_after_fork)


def new(api_key: str, base_url: str, timeout: int = 10) -> client.Client:
    """
    创建一个新的 SwanLab 运行时客户端。

    在共享锁内原子地检查、创建并赋值单例，覆盖网络认证构造期间，
    确保两个线程不能同时通过 exists() 检查并分别创建 Client。

    :return: 新创建的 Client 实例（所有权凭证）
    """
    global _default_client
    console.debug("Creating new SwanLab client.")
    with _client_lock:
        if _default_client is not None:
            raise RuntimeError("SwanLab client already exists. Call `reset` first.")
        _default_client = client.Client(api_key=api_key, base_url=base_url, timeout=timeout)
        return _default_client


def exists() -> bool:
    """检查当前的 SwanLab 运行时客户端是否已存在。"""
    global _default_client
    with _client_lock:
        return _default_client is not None


def reset(client: Optional[client.Client] = None):
    """
    重置/销毁当前的 SwanLab 运行时客户端，释放底层连接池。

    :param client: 调用方拥有的 Client 实例（所有权凭证）。
                          仅当全局单例仍是此实例时才执行销毁，避免误杀他人创建的实例。
                          None 时兼容旧代码（无所有权隔离），但不推荐。
    """
    global _default_client
    console.debug("Resetting SwanLab client.")
    with _client_lock:
        if _default_client is None:
            # 幂等：已被销毁，不报错
            return
        # 所有权隔离：仅当全局单例仍是调用方拥有的实例时才销毁
        if client is not None and _default_client is not client:
            console.debug(
                "Skipping client reset: global client instance has changed since creation. "
                "Another Core may have replaced it."
            )
            return
        _default_client.close()
        _default_client = None


def _get_client() -> client.Client:
    """获取当前的默认客户端。"""
    if _default_client is None:
        raise RuntimeError("SwanLab client is not initialized. Call `new` first.")
    return _default_client


def get(url: str, params: Optional[JSONDict] = None, retries: Optional[int] = None, log_error: bool = True):
    return _get_client().get(url, params=params, retries=retries, log_error=log_error)


def post(url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
    return _get_client().post(url, data=data, retries=retries, log_error=log_error)


def put(url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
    return _get_client().put(url, data=data, retries=retries, log_error=log_error)


def patch(url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
    return _get_client().patch(url, data=data, retries=retries, log_error=log_error)


def delete(url: str, retries: Optional[int] = None, log_error: bool = True):
    return _get_client().delete(url, retries=retries, log_error=log_error)


def username() -> str:
    return _get_client().username
