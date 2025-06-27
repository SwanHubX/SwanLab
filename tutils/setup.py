#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/27 14:57
@File: setup.py
@IDE: pycharm
@Description:
    存储、设置通用函数
"""
import os.path
from datetime import datetime, timedelta
from typing import Literal, Optional

import nanoid
import requests_mock

from swanlab.core_python import auth, create_client, reset_client, Client
from swanlab.data.store import reset_run_store, get_run_store, RunStore
from swanlab.package import get_host_web
from .config import TEMP_PATH

__all__ = ["UseMockRunState", "mock_login_info"]


class UseMockRunState:
    """
    使用上下文管理器模拟客户端以及一些其他的运行时状态
    主要用于测试环境中，模拟登录状态以及方便测试
    """

    def __init__(self):
        self.client: Optional[Client] = None
        self.store: Optional[RunStore] = None

    def __enter__(self) -> "UseMockRunState":
        login_info = mock_login_info()
        self.client = create_client(login_info)
        reset_run_store()
        self.store = get_run_store()
        # 创建运行目录结构，方便测试
        self.store.swanlog_dir = TEMP_PATH
        run_id = nanoid.generate("0123456789abcdefghijklmnopqrstuvwxyz", 21)
        self.store.run_id = run_id
        self.store.run_dir = os.path.join(TEMP_PATH, "run-" + run_id)
        os.mkdir(self.store.run_dir)
        os.mkdir(self.store.media_dir)
        os.mkdir(self.store.log_dir)
        os.mkdir(self.store.console_dir)
        os.mkdir(self.store.file_dir)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__del__()

    def __del__(self):
        reset_client()
        self.client = None
        reset_run_store()
        self.store = None


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
