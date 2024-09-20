#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/21 16:57
@File: pytest_log.py
@IDE: pycharm
@Description:
    测试swanlog类，只需测试其日志监听功能
"""
import os
import sys
import pytest
from swanlab.log import swanlog
from tutils import TEMP_PATH
from nanoid import generate
from freezegun import freeze_time


@pytest.fixture(scope="function", autouse=True)
def before_test_global_swanlog():
    """
    在测试之前清除全局swanlog的状态
    """
    try:
        swanlog.uninstall()
    except RuntimeError:
        pass


class TestSwanLogInstall:
    """
    目前在设计上不希望外界实例化SwanLog，所以不提供实例化测试
    """

    @staticmethod
    def create_console_dir():
        # 创建console文件夹
        console_dir = os.path.join(TEMP_PATH, str(generate()))
        os.mkdir(console_dir)
        return console_dir

    def test_global_install(self):
        swanlog.install()
        assert swanlog.installed is True
        with pytest.raises(RuntimeError) as e:
            swanlog.install()
        assert str(e.value) == "SwanLog has been installed"
        swanlog.uninstall()
        swanlog.install()
        assert swanlog.installed is True

    def test_write_after_uninstall(self):
        console_dir = self.create_console_dir()
        swanlog.install(console_dir)
        swanlog.uninstall()
        # 加一行防止其他问题
        print("\ntest write after uninstall")
        a = generate()
        print(a)
        b = generate()
        print(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        with open(os.path.join(console_dir, files[0]), "r") as f:
            content = f.readlines()
            assert len(content) == 0

    def test_write_to_file(self):
        console_dir = self.create_console_dir()
        swanlog.install(console_dir)
        # 加一行防止其他问题
        print("\ntest write to file")
        a = generate()
        print(a)
        b = generate()
        print(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        # 比较最后两行内容
        with open(os.path.join(console_dir, files[0]), "r") as f:
            content = f.readlines()
            assert content[-2] == a + "\n"
            assert content[-1] == b + "\n"

    def test_write_to_file_long_test(self):
        console_dir = self.create_console_dir()
        swanlog.install(console_dir)
        # 加一行防止其他问题
        print("\ntest write to file")
        a = generate(size=201)
        print(a)
        files = os.listdir(console_dir)
        with open(os.path.join(console_dir, files[0]), "r") as f:
            content = f.readlines()
            assert content[-1] == a[:200] + "\n"

    def test_write_logging_to_file(self):
        console_dir = self.create_console_dir()
        swanlog.install(console_dir, log_level="debug")
        # 加一行防止其他问题
        print("\ntest write to file")
        a = generate()
        swanlog.debug(a)
        b = generate()
        swanlog.info(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        with open(os.path.join(console_dir, files[0]), "r") as f:
            content = f.readlines()
            assert content[-2] == "swanlab: " + a + "\n"
            assert content[-1] == "swanlab: " + b + "\n"

    def test_can_write_logging(self):
        console_dir = self.create_console_dir()
        swanlog.install(console_dir)
        # 加一行防止其他问题
        print("\ntest write to file")
        a = generate()
        swanlog.debug(a)
        b = generate()
        swanlog.info(b)
        files = os.listdir(console_dir)
        assert len(files) == 1
        with open(os.path.join(console_dir, files[0]), "r") as f:
            content = f.readlines()
            assert content[-2] != "swanlab: " + a + "\n"
            assert content[-1] == "swanlab: " + b + "\n"

    def test_write_sharding(self, monkeypatch):
        """
        测试日志文件分片
        """
        console_dir = self.create_console_dir()
        with freeze_time('2020-10-06'):
            swanlog.install(console_dir)
            print("1234")
            assert os.path.exists(os.path.join(console_dir, "2020-10-06.log"))
        with freeze_time('2020-10-07'):
            p = os.path.join(console_dir, "2020-10-07.log")
            assert not os.path.exists(p)
            print("1234")
            assert os.path.exists(p)
