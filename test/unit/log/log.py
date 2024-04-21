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
from tutils import clear, SWANLAB_LOG_DIR
from nanoid import generate
import os


@pytest.fixture(scope="function", autouse=True)
def before_test_global_swanlog():
    """
    在测试之前清除全局swanlog的状态
    """
    try:
        swanlog.uninstall()
    except RuntimeError:
        pass
    # clear操作需要在uninstall之后
    clear()
    yield


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

    def test_write_after_uninstall(self, before_test_global_swanlog):
        console_dir = os.path.join(SWANLAB_LOG_DIR, 'console')
        os.mkdir(console_dir)
        swanlog.install(console_dir)
        swanlog.uninstall()
        # 加一行防止其他问题
        import sys
        # 开启标准输出流
        sys.stdout = sys.__stdout__
        print('\ntest write after uninstall')
        a = generate()
        print(a)
        b = generate()
        print(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        with open(os.path.join(console_dir, files[0]), 'r') as f:
            content = f.readlines()
            assert len(content) == 0
