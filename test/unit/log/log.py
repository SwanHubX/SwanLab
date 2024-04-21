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


class TestSwanLogInstall:
    def test_install_success(self):
        lg = SwanLog('tmp')
        lg.install()
        assert lg.installed is True
        lg.uninstall()

    def test_install_duplicate(self):
        lg = SwanLog('tmp')
        lg.install()
        with pytest.raises(RuntimeError) as e:
            lg.install()
        assert str(e.value) == 'SwanLog has been installed'
        lg.uninstall()

    def test_global_install(self, before_test_global_swanlog):
        swanlog.install()
        assert swanlog.installed is True
        with pytest.raises(RuntimeError) as e:
            swanlog.install()
        assert str(e.value) == 'SwanLog has been installed'

    def test_write_to_file(self, before_test_global_swanlog):
        console_dir = os.path.join(SWANLAB_LOG_DIR, 'console')
        os.mkdir(console_dir)
        swanlog.install(console_dir)
        # 加一行防止其他问题
        print('\ntest write to file')
        a = generate()
        print(a)
        b = generate()
        print(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        # 比较最后两行内容
        with open(os.path.join(console_dir, files[0]), 'r') as f:
            content = f.readlines()
            assert content[-2] == a + '\n'
            assert content[-1] == b + '\n'
        swanlog.uninstall()

    # def test_write_after_uninstall(self, before_test_global_swanlog):
    #     console_dir = os.path.join(SWANLAB_LOG_DIR, 'console')
    #     os.mkdir(console_dir)
    #     swanlog.install(console_dir)
    #     swanlog.uninstall()
    #     # 加一行防止其他问题
    #     print('\ntest write after uninstall')
    #     a = generate()
    #     print(a)
    #     b = generate()
    #     print(b)
    #     files = os.listdir(console_dir)
    #     assert len(files) == 0
