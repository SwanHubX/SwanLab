"""
@author: cunyue
@file: session.py
@time: 2025/9/9 15:10
@description: 创建会话
"""

import copy
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from swanlab.package import get_package_version

RETRY_HEADER = "X-Custom-Retry"
# 设置默认超时时间为 60s
DEFAULT_TIMEOUT = 60


class TimeoutHTTPAdapter(HTTPAdapter):
    """
    创建一个自定义的 HTTPAdapter，用于注入默认超时
    并可以通过请求 headers 中的 {RETRY_HEADER} 指定并临时修改重试次数
    """

    def __init__(self, *args, **kwargs):
        # 从 kwargs 中取出默认超时时间，如果没有则设为 None
        self.timeout = kwargs.pop("timeout", None)
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        # 如果 kwargs 中没有显式设置 timeout，则使用 self.timeout
        if "timeout" not in kwargs and self.timeout is not None:
            kwargs["timeout"] = self.timeout

        # 检查 headers 中是否有用户注入的重试次数
        _retry = request.headers.pop(RETRY_HEADER, None)
        if _retry is not None:
            _adapter = copy.copy(self)
            try:
                retry_count = int(_retry)
                if retry_count < 0:
                    raise ValueError("Retry count must be a non-negative integer.")
                _adapter.max_retries = self.max_retries.new(total=retry_count)
            except ValueError:
                raise ValueError(
                    f"Invalid retry count in {RETRY_HEADER}: '{_retry}'. Must be a non-negative integer."
                ) from None

            return _adapter.send(request, **kwargs)

        return super().send(request, **kwargs)


class SessionWithRetry(requests.Session):
    """
    自定义会话，用于自定义会话重试次数
    可以接受一个 retries 参数，并将其放在 headers 中的 {RETRY_HEADER} 字段中
    """

    def request(self, method, url, *args, **kwargs):
        retries = kwargs.pop('retries', None)

        # 将用户指定的重试次数注入到 headers 中
        if retries is not None:
            kwargs.setdefault('headers', {})[RETRY_HEADER] = str(retries)

        return super().request(method, url, *args, **kwargs)

    # ---------------------------------- 重写方法的函数签名，避免IDE警告 ----------------------------------

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


def create_session() -> SessionWithRetry:
    """
    创建一个带重试机制的会话
    :return: requests.Session
    """
    session = SessionWithRetry()
    retry = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET", "POST", "PUT", "DELETE", "PATCH"]),
    )
    adapter = TimeoutHTTPAdapter(max_retries=retry, timeout=DEFAULT_TIMEOUT)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers["swanlab-sdk"] = get_package_version()
    return session
