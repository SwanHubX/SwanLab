#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:30
@File: test_cli_utils.py
@IDE: pycharm
@Description:
    测试cli相关工具函数
"""
import pytest
from swanlab.cli.utils import UseTaskHttp
import tutils.setup as SU


def test_use_task_http_ok():
    with SU.UseMocker() as m:
        m.post("/test", text="mock")
        with SU.UseSetupHttp():
            with UseTaskHttp() as http:
                text = http.post("/test")
                assert text == "mock"


def test_use_task_http_abandon():
    with pytest.raises(SystemExit) as p:
        with SU.UseMocker() as m:
            m.post("/test", status_code=301, reason="Abandon")
            with SU.UseSetupHttp():
                with UseTaskHttp() as http:
                    http.post("/test")
    assert p.value.code == 3
