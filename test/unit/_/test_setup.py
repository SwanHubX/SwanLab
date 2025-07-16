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
    from swanlab.data.store import get_run_store

    with SU.UseMockRunState() as run_state:
        assert run_state.client is not None
        assert get_client() is not None
        assert run_state.store is not None
        run_store = get_run_store()
        run_store.run_dir = "1"
        assert run_state.store.run_dir == "1"

    assert get_run_store().run_dir is None, "Run store should be reset after UseMockRunState context manager exits"
    with pytest.raises(ValueError):
        get_client()
