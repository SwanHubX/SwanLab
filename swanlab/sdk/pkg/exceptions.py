"""
@author: cunyue
@file: exceptions.py
@time: 2026/3/7 21:36
@description: SwanLab 运行时客户端异常定义
"""

from typing import Union

from requests.exceptions import HTTPError


class ApiError(HTTPError):
    """
    SwanLab API 请求异常。
    封装了从后端解析出的业务错误码 (code) 和错误信息 (message)。
    """

    def __init__(self, response, code: Union[int, str], message: str, trace_id: str):
        self.response = response
        self.code = code
        self.message = message
        self.trace_id = trace_id
        self.request = response.request

        method = (self.request.method or "UNKNOWN").upper()

        # 构造友好的报错信息，控制台打印时一目了然
        error_str = f"API Request Failed: [{code}] {message} | Trace ID: {trace_id} | {method} {response.url}"
        super().__init__(error_str, response=response)
