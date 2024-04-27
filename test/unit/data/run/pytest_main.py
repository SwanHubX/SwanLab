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
import pytest


class TestSwanLabRun:

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
        assert str(e.value) == "SwanLabRun has been initialized"
        assert run.__str__() == get_run().__str__()
        _run = run.finish()
        assert swanlog.installed is False
