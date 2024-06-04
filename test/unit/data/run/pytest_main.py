#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_main.py
@IDE: pycharm
@Description:
    测试SwanLabRun主类
"""
import math

from swanlab.data.run.main import SwanLabRun, get_run, SwanLabRunState, swanlog
from swanlab import Image, Audio, Text
from nanoid import generate
from tutils import clear, TEMP_PATH
from PIL import Image as PILImage
import torch
import soundfile as sf
import numpy as np
import pytest
import random
import os


@pytest.fixture(scope="class", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    if get_run() is not None:
        get_run().finish()
    swanlog.disable_log()
    yield
    clear()
    swanlog.enable_log()


class TestSwanLabRunInit:

    def test_before_init(self):
        run = get_run()
        assert run is None
        assert SwanLabRun.get_state() == SwanLabRunState.NOT_STARTED

    def test_after_init(self):
        run = SwanLabRun(generate())
        assert swanlog.installed is False
        assert run is not None
        assert get_run().__str__() == run.__str__()
        _run = run.finish()
        assert get_run() is None
        assert _run.__str__() == run.__str__()
        assert SwanLabRunState.SUCCESS == _run.state

    def test_duplicate_init(self):
        run = SwanLabRun(generate())
        with pytest.raises(RuntimeError) as e:
            SwanLabRun(generate())
        assert swanlog.installed is False
        assert str(e.value) == "SwanLabRun has been initialized"
        assert run.__str__() == get_run().__str__()
        _run = run.finish()
        assert swanlog.installed is False


class TestSwanLabRunState:
    """
    测试SwanLabRun的状态变化
    """

    def test_not_started(self):
        assert SwanLabRun.get_state() == SwanLabRunState.NOT_STARTED

    def test_running(self):
        run = SwanLabRun()
        assert run.state == SwanLabRunState.RUNNING
        assert run.is_running is True

    def test_crashed(self):
        run = SwanLabRun()
        run.finish(SwanLabRunState.CRASHED, error="error")
        assert run.state == SwanLabRunState.CRASHED
        assert run.is_crashed is True

    def test_success(self):
        run = SwanLabRun()
        run.finish(SwanLabRunState.SUCCESS)
        assert run.state == SwanLabRunState.SUCCESS
        assert run.is_success is True


class TestSwanLabRunLog:
    """
    测试SwanLabRun的日志解析功能，不包含操作员
    """
    run = None

    @classmethod
    def setup_class(cls):
        cls.run = SwanLabRun()
        swanlog.disable_log()

    @classmethod
    def teardown_class(cls):
        cls.run.finish()
        swanlog.enable_log()

    # ---------------------------------- 解析log数字 ----------------------------------

    def test_log_number_ok(self):
        ll = self.run.log({"a": 1, "b": 0.1, "math.nan": math.nan, "math.inf": math.inf})
        assert len(ll) == 4
