#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/27 14:57
@File: setup.py
@IDE: pycharm
@Description:
    存储、设置通用函数
"""
import requests_mock
from swanlab.api import LoginInfo
from typing import Literal
from datetime import datetime, timedelta
import nanoid
from swanlab.package import get_host_api

__all__ = ["mock_login_info", "UseSetupHttp", "UseMocker"]


def mock_login_info(
        username=None,
        key=None,
        error_reason: Literal["OK", "Unauthorized", "Authorization Required", "Forbidden"] = "OK"
) -> LoginInfo:
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
    from swanlab.api.auth.login import login_request
    with requests_mock.Mocker() as m:
        if error_reason != "OK":
            if error_reason == "Authorization Required" or error_reason == "Unauthorized":
                status_code = 401
            elif error_reason == "Forbidden":
                status_code = 403
            else:
                status_code = 500
            m.post(f"{get_host_api()}/login/api_key", status_code=status_code, reason=error_reason)
        else:
            expired_at = datetime.now().isoformat()
            expired_at = (datetime.fromisoformat(expired_at) + timedelta(days=7)).isoformat() + 'Z'
            m.post(f"{get_host_api()}/login/api_key", json={
                "sid": nanoid.generate(),
                # 时间为当前时间加7天
                "expiredAt": expired_at,
                "userInfo": {
                    "username": username
                }
            }, status_code=200)
        resp = login_request(key)
        login_info = LoginInfo(resp, key)
    return login_info


class UseSetupHttp:
    """
    用于全局使用的http对象
    使用with关键字，自动登录，退出时自动重置http
    也可以使用del手动释放
    """

    def __init__(self):
        self.http = None

    def __enter__(self):
        from swanlab.api import create_http
        login_info = mock_login_info()
        self.http = create_http(login_info)
        return self.http

    def __exit__(self, exc_type, exc_val, exc_tb):
        from swanlab.api.http import reset_http
        reset_http()
        self.http = None

    def __del__(self):
        if self.http is not None:
            from swanlab.api.http import reset_http
            reset_http()
            self.http = None


class UseMocker(requests_mock.Mocker):
    """
    使用request_mock库进行mock测试，由于现在绝大部分请求都在get_host_api上，所以封装一层
    """

    def __init__(self, base_url: str = None):
        super().__init__()
        base_url = base_url or get_host_api()
        self.base_url = base_url

    def get(self, router, *args, **kwargs):
        return super().get(*(self.base_url + router, *args), **kwargs)

    def post(self, router, *args, **kwargs):
        return super().post(*(self.base_url + router, *args), **kwargs)

    def put(self, router, *args, **kwargs):
        return super().put(*(self.base_url + router, *args), **kwargs)

    def patch(self, router, *args, **kwargs):
        return super().patch(*(self.base_url + router, *args), **kwargs)

    def delete(self, router, *args, **kwargs):
        return super().delete(*(self.base_url + router, *args), **kwargs)
