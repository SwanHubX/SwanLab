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
    TODO 编写和媒体相关的测试脚本
    """

    @pytest.mark.parametrize("data", [
        random.randint(1, 100),
        -random.randint(1, 100),
        random.random(),
    ])
    def test_log_number_ok(self, data):
        """
        测试解析一个正常的数字
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": data})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert len(a.metric) == 3
        assert 'index' in a.metric
        assert 'create_time' in a.metric
        assert 'data' in a.metric
        assert a.metric['data'] == data
        assert a.step == 0
        assert a.epoch == 1
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_nan(self):
        """
        测试解析一个nan
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": float("nan")})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert a.metric is None
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error['data_class'] == "NaN"
        assert ac.error['excepted'] == ['float', 'int']
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_str(self):
        """
        测试解析其他字符串
        """
        """
        测试解析一个nan
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": 'abc'})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert a.metric is None
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error['data_class'] == "str"
        assert ac.error['excepted'] == ['float', 'int']
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_use_step(self):
        """
        测试解析一个数字，使用step
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": 1}, step=1)
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert len(a.metric) == 3
        assert 'index' in a.metric
        assert 'create_time' in a.metric
        assert 'data' in a.metric
        assert a.metric['data'] == 1
        assert a.step == 1
        assert a.epoch == 1
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_use_prefix(self):
        """
        测试解析一个数字，使用prefix
        """
        run = SwanLabRun()
        prefix_1 = generate()
        prefix_2 = generate() + '/' + generate()
        key1 = f"{prefix_1}/a"
        key2 = f"{prefix_2}/a"
        metric_dict = run.log({key1: 1, key2: 1})
        assert len(metric_dict) == 2
        for key in metric_dict:
            assert metric_dict[key] is not None
            a = metric_dict[key]
            ac = a.column_info
            # ---------------------------------- 指标信息 ----------------------------------
            assert len(a.metric) == 3
            assert 'index' in a.metric
            assert 'create_time' in a.metric
            assert 'data' in a.metric
            assert a.metric['data'] == 1
            assert a.step == 0
            assert a.epoch == 1
            # ---------------------------------- 列信息 ----------------------------------
            assert ac.data_type == "default"
            assert ac.error is None
            # 默认排在最前面
            assert ac.sort is None
            assert ac.config == {}
            assert ac.key == key
            assert ac.chart_type == "default"
            assert ac.namespace == key.split('/')[0]
            assert ac.reference == "step"
