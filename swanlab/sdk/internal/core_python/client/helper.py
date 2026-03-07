"""
@author: cunyue
@file: helper.py
@time: 2026/3/7 17:49
@description: SwanLab 运行时客户端辅助函数
"""

import json
from typing import Dict, List, Optional, Tuple, Union

import requests

from swanlab.sdk.pkg.helper import catch_and_return_none


def decode_response(resp: requests.Response) -> Union[Dict, List, str]:
    """
    解码响应，合并异常捕获
    """
    try:
        return resp.json()
    except (json.decoder.JSONDecodeError, requests.JSONDecodeError):
        return resp.text


@catch_and_return_none()
def decode_error_response(resp: requests.Response) -> Optional[Tuple[str, str]]:
    """
    尝试从错误响应的 JSON body 中解码业务错误码 (code) 和错误信息 (message)。
    如果响应为空、非 JSON 格式或缺少字段，装饰器会捕获异常并返回 None。

    :param resp: 错误响应
    :return: (code, message) 元组。如果解析失败则返回 None
    """
    # 如果响应体为空，提前结束
    if not resp.text.strip():
        return None

    data = resp.json()

    # 确保后端返回的是字典格式，并且包含了我们需要的键
    if isinstance(data, dict):
        # 即使后端没有严格同时返回 code 和 message，只要有其中之一也可以尽量提取
        # 提取不到的可以用原生的 status_code 和 reason 补位
        code = str(data.get("code", resp.status_code))
        message = str(data.get("message", resp.reason))
        return code, message

    return None
