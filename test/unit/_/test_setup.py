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
    from swanlab.core_python import get_client

    with SU.UseSetupHttp() as http:
        assert http is not None
        assert get_client() is not None
    with pytest.raises(ValueError):
        get_client()
