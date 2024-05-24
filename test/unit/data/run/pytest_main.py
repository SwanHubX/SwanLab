#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_main.py
@IDE: pycharm
@Description:
    测试SwanLabRun主类
"""
from swanlab.data.run.main import SwanLabRun, get_run, SwanLabRunState, swanlog
from nanoid import generate
from tutils import clear
import pytest
import random


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    if get_run() is not None:
        get_run().finish()
    yield
    clear()


class TestSwanLabRunInit:

    def test_before_init(self):
        run = get_run()
        assert run is None
        assert SwanLabRun.get_state() == SwanLabRunState.NOT_STARTED

    def test_after_init(self):
        run = SwanLabRun(generate())
        assert swanlog.installed is True
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
        assert swanlog.installed is True
        assert str(e.value) == "SwanLabRun has been initialized"
        assert run.__str__() == get_run().__str__()
        _run = run.finish()
        assert swanlog.installed is False


class TestSwanLabRunLog:
    """
    测试SwanLabRun的日志解析功能
    1. 输入为字典
    2. 可包含参数step
    3. 输入的字典的key必须为字符串
    4. 输出解析后的数据对象，包含额外的信息，step等信息
    5. 如果某一个解析失败，对应的key存在，但返回的数据为None
    """

    @pytest.mark.parametrize("data", [
        random.randint(1, 100),
        -random.randint(1, 100),
        random.random(),
    ])
    def test_log_number(self, data):
        """
        测试解析一个正常的数字
        """
        run = SwanLabRun()
        = run.log({"a": data})
