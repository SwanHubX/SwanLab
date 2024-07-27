#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/27 14:51
@File: test_cli_task.py
@IDE: pycharm
@Description:
    测试cli/task
"""
import pytest
from swanlab.cli.commands.task.utils import UseTaskHttp
import tutils.setup as SU
import requests_mock


def test_use_task_http_ok():
    with requests_mock.Mocker() as m:
        m.post("http://mock.swnalab.test/test", text="mock")
        with SU.UseSetupHttp():
            with UseTaskHttp() as http:
                resp = http.post("http://mock.swnalab.test/test")
                assert resp.text == "mock"


def test_use_task_http_abandon():
    with pytest.raises(SystemExit) as p:
        with SU.UseMocker() as m:
            m.post("/test", status_code=301, reason="Abandon")
            with SU.UseSetupHttp():
                with UseTaskHttp() as http:
                    http.post("/test")
    assert p.value.code == 3
