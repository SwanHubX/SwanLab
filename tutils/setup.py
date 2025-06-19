#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/27 14:57
@File: setup.py
@IDE: pycharm
@Description:
    存储、设置通用函数
"""
from datetime import datetime, timedelta
from typing import Literal

import nanoid
import requests_mock

from swanlab.core_python import auth
from swanlab.package import get_host_web

__all__ = ["UseSetupHttp", "mock_login_info"]


class UseSetupHttp:
    """
    用于全局使用的http对象，模拟登录，退出时重置http
    使用with关键字，自动登录，退出时自动重置http
    也可以使用del手动释放
    """

    def __init__(self):
        self.http = None

    def __enter__(self):
        from swanlab.core_python import create_client

        login_info = mock_login_info()
        self.http = create_client(login_info)
        return self.http

    def __exit__(self, exc_type, exc_val, exc_tb):
        from swanlab.core_python import reset_client

        reset_client()
        self.http = None

    def __del__(self):
        if self.http is not None:
            from swanlab.core_python import reset_client

            reset_client()
            self.http = None


def mock_login_info(
    username=None,
    key=None,
    error_reason: Literal["OK", "Unauthorized", "Authorization Required", "Forbidden"] = "OK",
) -> auth.LoginInfo:
    """
    生成一个虚假用户登录信息，主要用于本地mock，不能真实验证登录
    :param username: 需要使用的用户名，如果为None则随机生成
    :param key: 密钥，如果为None则随机生成
    :param error_reason: 错误原因,默认为OK，无错误
    :return: LoginInfo
    """
    if username is None:
        username = nanoid.generate()
    if key is None:
        key = nanoid.generate()
    from swanlab.package import get_host_api

    with requests_mock.Mocker() as m:
        api_host, web_host = get_host_api(), get_host_web()
        if error_reason != "OK":
            if error_reason == "Authorization Required" or error_reason == "Unauthorized":
                status_code = 401
            elif error_reason == "Forbidden":
                status_code = 403
            else:
                status_code = 500
            m.post(f"{api_host}/login/api_key", status_code=status_code, reason=error_reason)
        else:
            expired_at = datetime.now().isoformat()
            # 过期时间为当前时间加8天，主要是时区问题，所以不能7天以内
            expired_at = (datetime.fromisoformat(expired_at) + timedelta(days=8)).isoformat() + 'Z'
            m.post(
                f"{api_host}/login/api_key",
                json={"sid": nanoid.generate(), "expiredAt": expired_at, "userInfo": {"username": username}},
                status_code=200,
            )

        resp = auth.login_request(key, api_host)
        login_info = auth.LoginInfo(resp, key, api_host, web_host)
    return login_info
