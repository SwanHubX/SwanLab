#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:16
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试登录
"""
import os

import nanoid

from swanlab.env import SwanLabEnv
from swanlab.api.auth.login import login_by_key, terminal_login, code_login
from swanlab.error import ValidationError, APIKeyFormatError
from swanlab.package import is_login
from nanoid import generate
import tutils as T
import pytest


def get_password(prompt: str):
    # 如果是第一次登录，使用错误的key，会提示重新输入
    if "Paste" in prompt:
        return generate()
    else:
        return T.API_KEY


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_login_success():
    """
    测试登录成功
    """
    login_info = login_by_key(T.API_KEY, save=False)
    assert not login_info.is_fail
    assert login_info.api_key == T.API_KEY
    assert login_info.__str__() == "Login success"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_login_error_key():
    """
    测试登录失败, 错误的key
    """
    login_info = login_by_key("wrong-key", save=False)
    assert login_info.is_fail
    assert login_info.api_key is None
    assert login_info.__str__() == "Error api key"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_terminal_login(monkeypatch):
    """
    测试终端登录
    """
    monkeypatch.setattr("getpass.getpass", get_password)
    login_info = terminal_login(T.API_KEY)
    assert not login_info.is_fail
    assert login_info.api_key == T.API_KEY
    assert login_info.__str__() == "Login success"
    assert is_login()


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_code_login():
    """
    测试code登录
    """
    login_info = code_login(T.API_KEY)
    assert not login_info.is_fail
    assert login_info.api_key == T.API_KEY
    assert login_info.__str__() == "Login success"
    with pytest.raises(APIKeyFormatError):
        _ = code_login("wrong-key")


def test_api_key_format_error():
    """
    测试api_key格式错误: 类型、长度、非法字符
    """
    key = "123456"
    with pytest.raises(APIKeyFormatError) as e:
        _ = code_login(key)
    assert e.value.__str__() == "Api key length must be 21 characters long, yours was {}".format(len(key))
    key = nanoid.generate(size=20, alphabet="abcdefg") + '\x16'
    with pytest.raises(APIKeyFormatError) as e:
        _ = code_login(key)
    assert e.value.__str__() == "Invalid character in api key: '\\x16'"
    key = 123456
    with pytest.raises(APIKeyFormatError) as e:
        _ = code_login(key)  # type: ignore
    assert e.value.__str__() == "Api key must be a string"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_code_login_task_runtime():
    # task模式下不保存token
    os.environ[SwanLabEnv.RUNTIME.value] = 'task'
    login_info = code_login(T.API_KEY)
    assert not login_info.is_fail
    # 没有保存在本地
    del os.environ[SwanLabEnv.API_KEY.value]
    assert not is_login()
