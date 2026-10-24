"""
@author: cunyue
@file: exceptions.py
@time: 2026/3/7 21:36
@description: SwanLab 运行时异常定义
"""

from typing import Optional, Union

from requests.exceptions import HTTPError

__all__ = ["ApiError", "AuthenticationError", "DataStoreError"]


class ApiError(HTTPError):
    """
    SwanLab 请求异常。
    封装了从后端解析出的业务错误码 (code) 和错误信息 (message)；
    请求未发出或没有响应时（如协议校验失败）也允许直接构造，并通过 retryable 显式声明是否可重试。
    """

    def __init__(
        self,
        response=None,
        *,
        method: Optional[str] = None,
        trace_id: Optional[str] = None,
        code: Union[int, str] = "unknown code",
        message: str = "unknown error",
        retryable: Optional[bool] = None,
    ):
        self.response = response
        self.code = code
        self.message = message
        self.trace_id = trace_id
        self.request = getattr(response, "request", None)
        self._retryable = retryable
        # 构造友好的报错信息
        error_str = f"API Request Failed: [{code}] {message}"
        if trace_id is not None:
            error_str += f" | Trace ID: {trace_id}"
        if response is not None:
            error_str += f" | {method} {response.url}"
        super().__init__(error_str, response=response)

    @property
    def retryable(self) -> bool:
        """显式声明优先；否则 4xx 视为确定性失败，状态码缺失或 5xx 保留重试。"""
        if self._retryable is not None:
            return self._retryable
        status = getattr(self.response, "status_code", None)
        return status is None or status >= 500


class AuthenticationError(Exception):
    """
    SwanLab 认证失败异常。
    """

    pass


class DataStoreError(Exception):
    """
    DataStore 读写异常。
    """

    pass
