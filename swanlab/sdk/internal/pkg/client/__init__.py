"""
@author: cunyue
@file: __init__.py
@time: 2026/4/14 00:38
@description: SwanLab 客户端对象，与服务器进行交互
"""

from dataclasses import dataclass
from typing import Any, Optional

import requests

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.sdk.typings.pkg.client import JSONBody, JSONDict
from swanlab.sdk.typings.pkg.client.bootstrap import UserProfile

from .. import nrc
from . import session
from .utils import decode_response

__all__ = ["Client", "session", "ApiResponse", "decode_response", "verify_api_key"]


@dataclass
class ApiResponse:
    data: Any
    raw: requests.Response


class Client:
    """
    封装 SwanLab HTTP 请求的核心客户端对象。
    仅负责请求与重试拦截；鉴权凭据随 Authorization 头常驻，无会话刷新语义。
    """

    def __init__(self, api_key: str, base_url: str, timeout: int = 10):
        # 移除末尾的斜杠，防止 URL 拼接时出现双斜杠
        self._base_url = nrc.fmt(base_url) + "/api"
        self._profile: Optional[UserProfile] = None

        # 初始化时仅创建一次会话，以复用底层 TCP 连接池
        self._session: session.SessionWithRetry = session.create()
        # 挂载 ApiKey 常驻头，后续请求不再依赖任何临时会话
        self._session.headers["Authorization"] = f"ApiKey {api_key}"
        try:
            data = self.request("GET", "auth/verify", timeout=timeout, retries=0, log_error=False).data
            self._profile = _decode_profile(data)
        except ApiError as e:
            # 认证阶段失败必须释放连接池，避免资源泄漏
            self._session.close()
            if e.response.status_code in (401, 403):
                reason = e.message
                if reason == "unknown error":
                    reason = "invalid API key" if e.response.status_code == 401 else "access denied"
                raise AuthenticationError(f"Authentication failed: {reason}") from e
            if e.response.status_code == 404:
                raise RuntimeError(
                    "Backend version is too old for ApiKey authentication. "
                    "Please upgrade your SwanLab self-hosted deployment."
                ) from e
            raise
        except BaseException:
            self._session.close()
            raise

    @property
    def profile(self) -> UserProfile:
        assert self._profile is not None, "Client is not authenticated."
        return self._profile

    @property
    def username(self) -> str:
        return self.profile["username"]

    def close(self) -> None:
        """显式关闭底层会话，释放 TCP 连接池。"""
        self._session.close()

    # ---------------------------------- 实例 HTTP 方法 ----------------------------------
    def request(self, method: str, url: str, **kwargs):
        """基础请求方法封装"""
        full_url = f"{self._base_url}/{url.lstrip('/')}"
        resp = self._session.request(method, full_url, **kwargs)
        return ApiResponse(data=decode_response(resp), raw=resp)

    def get(self, url: str, params: Optional[JSONDict] = None, retries: Optional[int] = None, log_error: bool = True):
        return self.request("GET", url, params=params, retries=retries, log_error=log_error)

    def post(self, url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
        return self.request("POST", url, json=data, retries=retries, log_error=log_error)

    def put(self, url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
        return self.request("PUT", url, json=data, retries=retries, log_error=log_error)

    def patch(self, url: str, data: JSONBody = None, retries: Optional[int] = None, log_error: bool = True):
        return self.request("PATCH", url, json=data, retries=retries, log_error=log_error)

    def delete(self, url: str, retries: Optional[int] = None, log_error: bool = True):
        return self.request("DELETE", url, retries=retries, log_error=log_error)


def _decode_profile(data: Any) -> UserProfile:
    if not isinstance(data, dict):
        raise ValueError("Invalid authentication response: expected a user profile.")
    uid = data.get("uid")
    username = data.get("username")
    created_at = data.get("createdAt")
    avatar = data.get("avatar", "")
    name = data.get("name", "")
    verified = data.get("verified", False)
    if (
        type(uid) is not int
        or uid <= 0
        or not isinstance(username, str)
        or not username.strip()
        or not isinstance(created_at, str)
        or not isinstance(avatar, str)
        or not isinstance(name, str)
        or not isinstance(verified, bool)
    ):
        raise ValueError("Invalid authentication response: malformed user profile.")
    return UserProfile(uid=uid, username=username, createdAt=created_at, avatar=avatar, name=name, verified=verified)


def verify_api_key(base_url: str, api_key: str, timeout: int = 10) -> UserProfile:
    """
    一次性校验凭证并返回用户档案，失败透传异常。
    不驻留任何全局状态，供无副作用的临时探活使用。
    """
    client = Client(api_key=api_key, base_url=base_url, timeout=timeout)
    try:
        return client.profile
    finally:
        client.close()
