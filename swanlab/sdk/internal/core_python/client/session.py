"""
@author: cunyue
@file: session.py
@time: 2026/3/7 17:50
@description: SwanLab 运行时客户端会话辅助函数
具有默认重试次数和超时时间，也支持自定义重试次数和超时时间
"""

import contextvars
import copy
from typing import Optional

from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from swanlab.sdk.pkg.version import get_swanlab_version

__all__ = ["create", "TimeoutHTTPAdapter", "SessionWithRetry"]

VERSION_HEADER = "X-SwanLab-SDK-Version"
request_retries_ctx = contextvars.ContextVar("request_retries", default=None)


class TimeoutHTTPAdapter(HTTPAdapter):
    """
    支持默认超时的 HTTPAdapter。
    """

    def __init__(self, *args, **kwargs):
        self.timeout = kwargs.pop("timeout", None)
        super().__init__(*args, **kwargs)

    def send(self, request, *args, **kwargs):
        if self.timeout is not None:
            kwargs.setdefault("timeout", self.timeout)

        # 2. 直接从上下文中读取 retries，无需触碰 request.headers
        retries = request_retries_ctx.get()

        if retries is not None:
            if not isinstance(retries, int) or retries < 0:
                raise ValueError(f"Invalid retry count: '{retries}'. Must be a non-negative integer.")

            adapter = copy.copy(self)
            adapter.max_retries = self.max_retries.new(total=retries)
            return super(TimeoutHTTPAdapter, adapter).send(request, **kwargs)

        return super().send(request, *args, **kwargs)


class SessionWithRetry(Session):
    """
    支持在请求级别自定义重试次数的 Session。
    通过拦截 retries 参数并将其转化为隐式 Header 传递给 Adapter。
    """

    def request(self, method, url, *args, **kwargs):
        retries = kwargs.pop("retries", None)

        if retries is not None:
            # 3. 将自定义参数放入上下文，并获取 token 以便后续清理
            token = request_retries_ctx.set(retries)
            try:
                return super().request(method, url, *args, **kwargs)
            finally:
                # 4. 请求结束后，务必清理上下文，避免影响复用该线程的其他请求
                request_retries_ctx.reset(token)
        else:
            return super().request(method, url, *args, **kwargs)

    # ---------------------------------- 类型提示占位符，保留以保证 IDE 友好 ----------------------------------

    def get(self, url, params=None, retries: Optional[int] = None, **kwargs):
        return self.request("GET", url, params=params, retries=retries, **kwargs)

    def options(self, url, retries: Optional[int] = None, **kwargs):
        return self.request("OPTIONS", url, retries=retries, **kwargs)

    def head(self, url, retries: Optional[int] = None, **kwargs):
        return self.request("HEAD", url, retries=retries, **kwargs)

    def post(self, url, data=None, json=None, retries: Optional[int] = None, **kwargs):
        return self.request("POST", url, data=data, json=json, retries=retries, **kwargs)

    def put(self, url, data=None, retries: Optional[int] = None, **kwargs):
        return self.request("PUT", url, data=data, retries=retries, **kwargs)

    def patch(self, url, data=None, retries: Optional[int] = None, **kwargs):
        return self.request("PATCH", url, data=data, retries=retries, **kwargs)

    def delete(self, url, retries: Optional[int] = None, **kwargs):
        return self.request("DELETE", url, retries=retries, **kwargs)


def create(timeout: int = 60) -> SessionWithRetry:
    """
    创建一个挂载了超时和重试机制的会话实例。
    """
    session = SessionWithRetry()

    retry_strategy = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET", "POST", "PUT", "DELETE", "PATCH"]),
    )

    adapter = TimeoutHTTPAdapter(max_retries=retry_strategy, timeout=timeout)

    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers[VERSION_HEADER] = get_swanlab_version()
    return session
