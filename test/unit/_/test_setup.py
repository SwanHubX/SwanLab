#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/27 15:07
@File: setup.py
@IDE: pycharm
@Description:
    测试tutils/setup.py
"""
import pytest
import tutils.setup as SU


def test_mock_login_info():
    login_info = SU.mock_login_info()
    assert login_info.is_fail is False
    login_info = SU.mock_login_info(error_reason="Unauthorized")
    assert login_info.is_fail is True
    login_info = SU.mock_login_info(error_reason="Authorization Required")
    assert login_info.is_fail is True
    login_info = SU.mock_login_info(error_reason="Forbidden")
    assert login_info.is_fail is True
    login_info = SU.mock_login_info(error_reason="OK")
    assert login_info.is_fail is False


def test_use_setup_http():
    from swanlab.api import get_http
    with SU.UseSetupHttp() as http:
        assert http is not None
        assert get_http() is not None
    with pytest.raises(ValueError):
        get_http()


def test_use_mocker():
    with SU.UseMocker() as m:
        m.post("/tmp", text="mock")
        import requests
        from swanlab.package import get_host_api
        resp = requests.post(get_host_api() + "/tmp")
        assert resp.text == "mock"
