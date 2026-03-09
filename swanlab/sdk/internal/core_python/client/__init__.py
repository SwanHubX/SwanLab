"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:04
@description: SwanLab 运行时客户端，用于与 SwanLab API 进行交互
"""

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional, Union

import requests

from swanlab.sdk.internal.core_python.api.bootstrap import login_by_api_key
from swanlab.sdk.internal.pkg.netrc import remove_host_suffix
from swanlab.sdk.internal.pkg.scope import set_context
from swanlab.sdk.typings.core_python.api.bootstrap import LoginResponse
from swanlab.sdk.utils.version import get_swanlab_version

from ...pkg import console
from . import session
from .helper import decode_response

__all__ = [
    "Client",
    "session",
    "exists",
    "reset",
    "new",
    "get",
    "post",
    "put",
    "patch",
    "delete",
]


@dataclass
class ApiResponse:
    data: Union[dict, list, str]
    raw: requests.Response


class Client:
    """
    封装 SwanLab HTTP 请求的核心客户端对象。
    仅负责请求、重试拦截及鉴权生命周期管理。
    """

    # 提前刷新的缓冲时间（7天）
    REFRESH_TIME = 60 * 60 * 24 * 7

    def __init__(self, api_key: str, base_url: str, timeout: int = 10):
        self._api_key = api_key
        self._version = get_swanlab_version()
        # 移除末尾的斜杠，防止 URL 拼接时出现双斜杠
        self._base_url = remove_host_suffix(base_url, "/api") + "/api"
        self._expired_at: Optional[datetime] = None

        # 初始化时仅创建一次会话，以复用底层 TCP 连接池
        self._session: session.SessionWithRetry = session.create()
        # 立即进行首次鉴权并挂载凭证
        login_resp = self._refresh_auth(timeout=timeout, warning=False)
        # 写入登录响应到上下文，由调用者判断是否需要使用
        set_context("login_resp", login_resp)

    def _refresh_auth(self, timeout: int = 10, warning: bool = True) -> Optional[LoginResponse]:
        """
        刷新鉴权信息。
        直接更新当前会话的 Cookie，保留底层连接池以提升性能。
        """
        console.debug("Refreshing authentication token...")
        login_resp = login_by_api_key(self._base_url, self._api_key, timeout=timeout)
        if not login_resp:
            if warning:
                console.warning(
                    "Failed to refresh authentication token,",
                    "swanlab may not work properly, please check your API key.",
                )
            return None
        # 只需更新当前 session 的 cookie 即可
        self._session.cookies.update({"sid": login_resp["sid"]})
        self._expired_at = datetime.strptime(login_resp["expiredAt"], "%Y-%m-%dT%H:%M:%S.%fZ").replace(
            tzinfo=timezone.utc
        )
        return login_resp

    def _before_request(self):
        """请求前置检查。距过期时间不足安全缓冲期时，触发刷新。"""
        if self._expired_at is None:
            return self._refresh_auth()

        if (self._expired_at - datetime.now(timezone.utc)).total_seconds() <= self.REFRESH_TIME:
            console.debug("Session is about to expire. Triggering token refresh.")
            self._refresh_auth()
        return None

    # ---------------------------------- 实例 HTTP 方法 ----------------------------------
    def request(self, method: str, url: str, **kwargs):
        """基础请求方法封装"""
        full_url = f"{self._base_url}/{url.lstrip('/')}"
        self._before_request()

        resp = self._session.request(method, full_url, **kwargs)
        return ApiResponse(data=decode_response(resp), raw=resp)

    def get(self, url: str, params: Optional[dict] = None, retries: Optional[int] = None):
        return self.request("GET", url, params=params, retries=retries)

    def post(self, url: str, data: Optional[Union[dict, list]] = None, retries: Optional[int] = None):
        return self.request("POST", url, json=data, retries=retries)

    def put(self, url: str, data: Optional[dict] = None, retries: Optional[int] = None):
        return self.request("PUT", url, json=data, retries=retries)

    def patch(self, url: str, data: Optional[dict] = None, retries: Optional[int] = None):
        return self.request("PATCH", url, json=data, retries=retries)

    def delete(self, url: str, retries: Optional[int] = None):
        return self.request("DELETE", url, retries=retries)


# ==============================================================================
# 模块级全局状态与代理快捷函数
# ==============================================================================

_default_client: Optional[Client] = None


def new(api_key: str, base_url: str, timeout: int = 10) -> Client:
    """创建一个新的 SwanLab 运行时客户端。"""
    global _default_client
    console.debug("Creating new SwanLab client.")
    if _default_client is not None:
        raise RuntimeError("SwanLab client already exists. Call `reset` first.")
    _default_client = Client(api_key=api_key, base_url=base_url, timeout=timeout)
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


def _get_client() -> Client:
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
