"""
@author: cunyue
@file: bootstrap.py
@time: 2026/3/7 18:30
@description: 这部分接口不经过Client调用，因为会在Client初始化之前调用
"""

from swanlab.sdk.internal.core_python.client import session
from swanlab.sdk.typings.core_python.api.bootstrap import LoginResponse


def login_by_api_key(base_url: str, api_key: str, timeout: int = 20) -> LoginResponse:
    """
    用户登录，请求后端接口完成验证
    """
    with session.create() as s:
        resp = s.post(url=f"{base_url}/login/api_key", headers={"authorization": api_key}, timeout=timeout)
    return resp.json()
