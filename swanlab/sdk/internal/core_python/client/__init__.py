"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:04
@description: SwanLab 运行时客户端，用于与 SwanLab API 进行交互
"""

from typing import Optional, Union

from swanlab.sdk.internal.pkg import client, console

__all__ = ["exists", "reset", "new", "get", "post", "put", "patch", "delete"]

# ==============================================================================
# 模块级全局状态与代理快捷函数
# ==============================================================================

_default_client: Optional[client.Client] = None


def new(api_key: str, base_url: str, timeout: int = 10) -> client.Client:
    """创建一个新的 SwanLab 运行时客户端。"""
    global _default_client
    console.debug("Creating new SwanLab client.")
    if _default_client is not None:
        raise RuntimeError("SwanLab client already exists. Call `reset` first.")
    _default_client = client.Client(api_key=api_key, base_url=base_url, timeout=timeout)
    return _default_client


def exists() -> bool:
    """检查当前的 SwanLab 运行时客户端是否已存在。"""
    global _default_client
    return _default_client is not None


def reset():
    """重置/销毁当前的 SwanLab 运行时客户端。"""
    global _default_client
    console.debug("Resetting SwanLab client.")
    if _default_client is None:
        raise RuntimeError("SwanLab client is not initialized. Call `new` first.")
    _default_client = None


def _get_client() -> client.Client:
    """获取当前的默认客户端。"""
    if _default_client is None:
        raise RuntimeError("SwanLab client is not initialized. Call `new` first.")
    return _default_client


def get(url: str, params: Optional[dict] = None, retries: Optional[int] = None):
    return _get_client().get(url, params=params, retries=retries)


def post(url: str, data: Optional[Union[dict, list]] = None, retries: Optional[int] = None):
    return _get_client().post(url, data=data, retries=retries)


def put(url: str, data: Optional[dict] = None, retries: Optional[int] = None):
    return _get_client().put(url, data=data, retries=retries)


def patch(url: str, data: Optional[dict] = None, retries: Optional[int] = None):
    return _get_client().patch(url, data=data, retries=retries)


def delete(url: str, retries: Optional[int] = None):
    return _get_client().delete(url, retries=retries)
