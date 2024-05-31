#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:16
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试登录
"""
from swanlab.api.auth.login import login_by_key, terminal_login, code_login
from swanlab.error import ValidationError
from swanlab.package import is_login
import tutils as T
from swanlab.env import reset_env, HOME
import pytest
import os


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    reset_env()
    yield
    reset_env()
    if HOME in os.environ:
        del os.environ[HOME]


def test_login_success():
    """
    测试登录成功
    """
    login_info = login_by_key(T.KEY, save=False)
    assert not login_info.is_fail
    assert login_info.api_key == T.KEY
    assert login_info.__str__() == "Login success"


def test_login_error_key():
    """
    测试登录失败, 错误的key
    """
    login_info = login_by_key("wrong-key", save=False)
    assert login_info.is_fail
    assert login_info.api_key is None
    assert login_info.__str__() == "Error api key"


def test_terminal_login(monkeypatch):
    """
    测试终端登录
    """
    os.environ[HOME] = T.TEMP_PATH
    monkeypatch.setattr("getpass.getpass", T.get_password)
    login_info = terminal_login(T.KEY)
    assert not login_info.is_fail
    assert login_info.api_key == T.KEY
    assert login_info.__str__() == "Login success"
    assert is_login()


def test_code_login():
    """
    测试code登录
    """
    login_info = code_login(T.KEY)
    assert not login_info.is_fail
    assert login_info.api_key == T.KEY
    assert login_info.__str__() == "Login success"
    with pytest.raises(ValidationError):
        _ = code_login("wrong-key")
