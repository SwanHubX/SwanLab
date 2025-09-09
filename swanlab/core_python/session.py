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
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers["swanlab-sdk"] = get_package_version()
    return session
