"""
@author: cunyue
@file: bootstrap.py
@time: 2026/3/7 18:30
@description: 这部分接口不经过Client调用，因为会在Client初始化之前调用
"""

from typing import Optional

from requests.exceptions import RequestException

from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

from . import session
from .session import STARTUP_RETRIES


def login_by_api_key(
    base_url: str,
    api_key: str,
    timeout: int = 20,
    retries: Optional[int] = STARTUP_RETRIES,
) -> LoginResponse:
    """
    用户登录，请求后端接口完成验证。
    网络错误、HTTP 错误和响应解析错误会向上抛出；Client 首次认证会统一包装为
    AuthenticationError，并将原始异常保留在 __cause__ 中。

    :param retries: 覆盖会话默认重试次数。默认使用 STARTUP_RETRIES（启动期间有限重试），
        运行期刷新传 None 以使用会话默认重试。
    """
    with session.create(timeout=timeout, default_retry=retries) as s:
        resp = s.post(url=f"{base_url}/login/api_key", headers={"authorization": api_key}, timeout=timeout)
    return resp.json()


@safe.decorator(RequestException, message=None)
def verify_api_key(base_url: str, api_key: str, timeout: int = 20) -> Optional[LoginResponse]:
    """
    验证 API Key，失败时返回 None。
    仅用于 CLI 校验场景，不参与 init 登录降级。
    """
    return login_by_api_key(base_url, api_key, timeout=timeout)
