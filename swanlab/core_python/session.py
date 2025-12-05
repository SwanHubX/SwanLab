"""
@author: cunyue
@file: session.py
@time: 2025/9/9 15:10
@description: 创建会话
"""

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from swanlab.package import get_package_version

# 设置默认超时时间为 60s
DEFAULT_TIMEOUT = 60


# 创建一个自定义的 HTTPAdapter，用于注入默认超时
class TimeoutHTTPAdapter(HTTPAdapter):
    def __init__(self, *args, **kwargs):
        # 从 kwargs 中取出默认超时时间，如果没有则设为 None
        self.timeout = kwargs.pop("timeout", None)
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        # 如果 kwargs 中没有显式设置 timeout，则使用 self.timeout
        if "timeout" not in kwargs and self.timeout is not None:
            kwargs["timeout"] = self.timeout

        return super().send(request, **kwargs)


def create_session() -> requests.Session:
    """
    创建一个带重试机制的会话
    :return: requests.Session
    """
    session = requests.Session()
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
