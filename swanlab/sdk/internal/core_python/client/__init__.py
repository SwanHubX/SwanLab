"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:04
@description: SwanLab 运行时客户端，用于与 SwanLab API 进行交互
"""

from datetime import datetime, timezone
from typing import Optional, Union

import requests

from swanlab.sdk.internal.core_python.api.auth import login_by_api_key
from swanlab.sdk.pkg import console
from swanlab.sdk.pkg.version import get_swanlab_version

from . import session
from .helper import decode_response

__all__ = [
    "Client",
    "session",
    "init_client",
    "get_client",
    "reset_client",
    "get",
    "post",
    "put",
    "patch",
    "delete",
]


class Client:
    """
    封装 SwanLab HTTP 请求的核心客户端对象。
    仅负责请求、重试拦截及鉴权生命周期管理。
    """

    # 提前刷新的缓冲时间（原逻辑为7天）
    REFRESH_TIME = 60 * 60 * 24 * 7

    def __init__(self, api_key: str, base_url: str):
        self._api_key = api_key
        self._version = get_swanlab_version()
        # 移除末尾的斜杠，防止 URL 拼接时出现双斜杠
        self._base_url = base_url.rstrip("/")
        self._expired_at: Optional[datetime] = None

        # 1. 初始化时仅创建一次会话，以复用底层 TCP 连接池
        self._session: session.SessionWithRetry = session.create()
        self._setup_interceptor()

        # 2. 立即进行首次鉴权并挂载凭证
        self._refresh_auth()

    def _setup_interceptor(self):
        """挂载全局响应拦截器"""

        def response_interceptor(response: requests.Response):
            method = (response.request.method or "UNKNOWN").upper()
            console.debug(f"HTTP Request: {method} {response.url} | Response Status: {response.status_code}")
            if response.status_code // 100 != 2:
                # TODO 解析错误信息
                raise Exception(
                    response,
                    f"Trace id: {response.headers.get('traceid')} {method} {response.url}",
                    f"{response.status_code} {response.reason}",
                )

        self._session.hooks["response"] = response_interceptor

    def _refresh_auth(self):
        """
        刷新鉴权信息。
        直接更新当前会话的 Cookie，保留底层连接池以提升性能。
        """
        console.debug("Refreshing authentication token...")

        login_resp = login_by_api_key(self._base_url, self._api_key)

        # 核心修改：只需更新当前 session 的 cookie 即可
        self._session.cookies.update({"sid": login_resp["sid"]})

    def _before_request(self):
        """请求前置检查。距过期时间不足安全缓冲期时，触发刷新。"""
        if self._expired_at is None:
            self._refresh_auth()
            return

        if (self._expired_at - datetime.now(timezone.utc)).total_seconds() <= self.REFRESH_TIME:
            console.debug("Session is about to expire. Triggering token refresh.")
            self._refresh_auth()

    # ---------------------------------- 实例 HTTP 方法 ----------------------------------
    def request(self, method: str, url: str, **kwargs):
        """基础请求方法封装"""
        full_url = f"{self._base_url}/{url.lstrip('/')}"
        self._before_request()

        resp = self._session.request(method, full_url, **kwargs)
        return decode_response(resp), resp

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


def init_client(api_key: str, base_url: str) -> Client:
    """初始化全局默认客户端。"""
    global _default_client
    _default_client = Client(api_key=api_key, base_url=base_url)
    return _default_client


def get_client() -> Client:
    """获取当前的默认客户端。"""
    if _default_client is None:
        raise RuntimeError("SwanLab client is not initialized. Call `init_client` first.")
    return _default_client


def reset_client():
    """重置/销毁全局客户端。"""
    global _default_client
    _default_client = None


def get(url: str, params: Optional[dict] = None, retries: Optional[int] = None):
    return get_client().get(url, params=params, retries=retries)


def post(url: str, data: Optional[Union[dict, list]] = None, retries: Optional[int] = None):
    return get_client().post(url, data=data, retries=retries)


def put(url: str, data: Optional[dict] = None, retries: Optional[int] = None):
    return get_client().put(url, data=data, retries=retries)


def patch(url: str, data: Optional[dict] = None, retries: Optional[int] = None):
    return get_client().patch(url, data=data, retries=retries)


def delete(url: str, retries: Optional[int] = None):
    return get_client().delete(url, retries=retries)
