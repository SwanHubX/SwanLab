"""
@author: cunyue
@file: bootstrap.py
@time: 2026/3/7 18:30
@description: 这部分接口不经过Client调用，因为会在Client初始化之前调用
"""

from typing import Optional

from requests.exceptions import RequestException

from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

from . import session

# 登录握手的重试次数。启动路径需要快速失败，避免长时间阻塞 init()。
LOGIN_RETRIES = 1


def login_by_api_key(base_url: str, api_key: str, timeout: int = 20) -> LoginResponse:
    """
    用户登录，请求后端接口完成验证。
    网络错误、HTTP 错误和响应解析错误会向上抛出，由调用方区分认证失败与服务不可用。
    """
    with session.create(timeout=timeout, default_retry=LOGIN_RETRIES) as s:
        resp = s.post(url=f"{base_url}/login/api_key", headers={"authorization": api_key}, timeout=timeout)
    return resp.json()


def verify_api_key(base_url: str, api_key: str, timeout: int = 20) -> Optional[LoginResponse]:
    """
    验证 API Key，失败时返回 None。
    仅用于 CLI 校验场景，不参与 init 登录降级。
    """
    try:
        return login_by_api_key(base_url, api_key, timeout=timeout)
    except RequestException:
        return None
