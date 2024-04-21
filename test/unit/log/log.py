#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/21 16:57
@File: log.py
@IDE: pycharm
@Description:
    测试swanlog类
"""
import pytest
from swanlab.log import SwanLog, swanlog


class TestSwanLogInstall:
    def test_install_success(self):
        lg = SwanLog('tmp')
        lg.install()
        assert lg.installed is True

    def test_install_duplicate(self):
        lg = SwanLog('tmp')
        lg.install()
        with pytest.raises(RuntimeError) as e:
            lg.install()
        assert str(e.value) == 'SwanLog has been installed'

    def test_global_install(self):
        swanlog.install()
        assert swanlog.installed is True
        with pytest.raises(RuntimeError) as e:
            swanlog.install()
        assert str(e.value) == 'SwanLog has been installed'
