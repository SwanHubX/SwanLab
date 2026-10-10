"""
@author: cunyue
@file: helper.py
@time: 2026/3/7 17:49
@description: SwanLab 运行时客户端辅助函数
"""

import json
from typing import Any, Dict, List, Optional, Tuple, Union

import requests

from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.pkg.client.bootstrap import UserProfile


def decode_response(resp: requests.Response) -> Union[Dict, List, str]:
    """
    解码响应，合并异常捕获
    """
    try:
        return resp.json()
    except (json.decoder.JSONDecodeError, requests.JSONDecodeError):
        return resp.text


@safe.decorator(message=None)
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


def decode_profile(data: Any) -> UserProfile:
    """解码认证返回的用户身份信息。"""
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
