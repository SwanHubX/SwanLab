"""
@author: cunyue
@file: __init__.py
@time: 2026/4/14 00:38
@description: SwanLab 客户端对象，与服务器进行交互
"""

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Optional

import requests
from requests.exceptions import InvalidSchema, InvalidURL, MissingSchema

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.sdk.typings.pkg.client import JSONBody, JSONDict
from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

from .. import console, helper, nrc, scope
from . import session
from .bootstrap import login_by_api_key, verify_api_key
from .session import STARTUP_RETRIES
from .utils import decode_response

__all__ = [
    "Client",
    "session",
    "ApiResponse",
    "decode_response",
    "verify_api_key",
    "is_transient_login_error",
    "is_login_auth_error",
    "STARTUP_RETRIES",
]


def is_transient_login_error(exc: BaseException) -> bool:
    """
    判断登录异常是否属于服务不可用类错误。

    Client 首次认证会把底层异常统一包装为 ``AuthenticationError``，原始异常保留在 ``__cause__`` 中， 因此这里只解包一层即可判断真实原因。
    """
    # 这里只解包这一层，避免递归遍历异常链
    if isinstance(exc, AuthenticationError) and exc.__cause__ is not None:
        exc = exc.__cause__
    if isinstance(exc, ApiError):
        status = getattr(exc.response, "status_code", None)
        return status is None or status in (408, 429) or status >= 500
    # URL / 协议配置错误重试和降级都无法解决
    if isinstance(exc, (InvalidURL, MissingSchema, InvalidSchema)):
        return False
    return isinstance(exc, requests.RequestException)


def is_login_auth_error(exc: BaseException) -> bool:
    """
    判断异常是否由 API Key 认证失败（401/403）导致；这类错误可以重新输入 Key 解决。
    与 ``is_transient_login_error`` 一样，只解包 Client 包装的一层 ``__cause__``。
    """
    if isinstance(exc, AuthenticationError) and exc.__cause__ is not None:
        exc = exc.__cause__
    return isinstance(exc, ApiError) and getattr(exc.response, "status_code", None) in (401, 403)


@dataclass
class ApiResponse:
    data: Any
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
        self._version = helper.get_swanlab_version()
        # 移除末尾的斜杠，防止 URL 拼接时出现双斜杠
        self._base_url = nrc.fmt(base_url) + "/api"
        self._expired_at: Optional[datetime] = None
        self._username: Optional[str] = None

        # 初始化时仅创建一次会话，以复用底层 TCP 连接池
        self._session: session.SessionWithRetry = session.create()
        # 立即进行首次鉴权并挂载凭证；首次认证失败统一包装为 AuthenticationError，
        # 原始异常保留在 __cause__，供 init 判断是否可降级
        login_resp = self._refresh_auth(timeout=timeout, warning=False, strict=True)
        # 写入登录响应到上下文，由调用者判断是否需要使用
        scope.set_context("login_resp", login_resp)

    def _refresh_auth(self, timeout: int = 10, warning: bool = True, strict: bool = False) -> Optional[LoginResponse]:
        """
        刷新鉴权信息。
        直接更新当前会话的 Cookie，保留底层连接池以提升性能。

        :param strict: 首次认证时为 True，异常统一包装为 AuthenticationError 并向上抛出；
            运行中刷新为 False，失败仅告警并保留会话默认重试。
        """
        console.debug("Refreshing authentication token...")
        retries = STARTUP_RETRIES if strict else None
        try:
            login_resp = login_by_api_key(self._base_url, self._api_key, timeout=timeout, retries=retries)
            try:
                sid = login_resp["sid"]
                expired_at = login_resp["expiredAt"]
                username = login_resp["userInfo"]["username"]
                parsed_expired_at = datetime.strptime(expired_at, "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)
            except (KeyError, TypeError, ValueError) as e:
                # 200 但响应体为空/缺字段，按服务端异常处理，不误判为 API Key 错误
                raise ApiError(
                    None,
                    method="POST",
                    code="invalid_login_response",
                    message=f"Invalid login response: {e}",
                    retryable=True,
                ) from e
            # 只需更新当前 session 的 cookie 即可
            self._session.cookies.update({"sid": sid})
            self._expired_at = parsed_expired_at
            self._username = username
            return login_resp
        except ApiError as e:
            if strict:
                # 错误在发生处打印调用链；原始异常保留在 __cause__ 供 init 判断是否可降级
                console.trace("Failed to login")
                raise self._wrap_auth_error(e) from e
            self._warn_refresh_failure(warning)
            return None
        except Exception as e:
            if strict:
                console.trace("Failed to login")
                raise self._wrap_auth_error(e) from e
            self._warn_refresh_failure(warning)
            return None

    @staticmethod
    def _wrap_auth_error(exc: BaseException) -> AuthenticationError:
        """把首次认证的底层异常包装为公开的 AuthenticationError。"""
        status = getattr(getattr(exc, "response", None), "status_code", None)
        if status in (401, 403):
            return AuthenticationError(
                "Failed to initialize the SwanLab client. Please check if the provided API key is correct."
            )
        return AuthenticationError(f"Failed to initialize the SwanLab client: {exc}")

    @staticmethod
    def _warn_refresh_failure(warning: bool) -> None:
        if not warning:
            return
        console.warning(
            "Failed to refresh authentication token,",
            "swanlab may not work properly, please check your API key.",
        )

    def _before_request(self):
        """请求前置检查。距过期时间不足安全缓冲期时，触发刷新。"""
        if self._expired_at is None:
            return self._refresh_auth()

        if (self._expired_at - datetime.now(timezone.utc)).total_seconds() <= self.REFRESH_TIME:
            console.debug("Session is about to expire. Triggering token refresh.")
            self._refresh_auth()
        return None

    @property
    def username(self) -> str:
        assert self._username is not None, "Username is not set. Ensure that the client is authenticated."
        return self._username

    # ---------------------------------- 实例 HTTP 方法 ----------------------------------
    def request(self, method: str, url: str, **kwargs):
        """基础请求方法封装"""
        full_url = f"{self._base_url}/{url.lstrip('/')}"
        self._before_request()

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
