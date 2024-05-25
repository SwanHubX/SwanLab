#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/24 20:25
@File: pytest_env.py
@IDE: pycharm
@Description:
    测试swanlab/env.py中的环境变量和工具函数
"""
from swanlab.env import (
    SwanLabMode,
    get_mode,
    assert_exist,
    reset_env,
    is_strict_mode,
    MODE
)
from nanoid import generate
import pytest
import os


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    reset_env()
    if MODE in os.environ:
        del os.environ[MODE]
    yield
    reset_env()
    if MODE in os.environ:
        del os.environ[MODE]


def use_strict_mode():
    """
    使用严格模式
    """
    get_mode({MODE: SwanLabMode.CLOUD.value})


def use_not_strict_mode():
    """
    使用非严格模式
    """
    get_mode({MODE: SwanLabMode.DISABLED.value})


class TestStrictMode:

    def test_strict_mode(self):
        """
        测试默认严格模式
        """
        assert is_strict_mode() is True

    def test_strict_mode_cloud(self):
        """
        测试严格模式
        """
        use_strict_mode()
        assert is_strict_mode() is True

    def test_strict_mode_disabled(self):
        """
        测试非严格模式
        """
        use_not_strict_mode()
        assert is_strict_mode() is False

    def test_strict_mode_local(self):
        """
        测试本地模式
        """
        use_strict_mode()
        assert is_strict_mode() is True

    def test_strict_mode_error(self):
        """
        测试重定向模式
        """
        with pytest.raises(ValueError):
            get_mode({MODE: generate()})


class TestAssertExist:

    def test_exist_ok(self):
        """
        测试严格模式下，文件路径存在
        """
        assert assert_exist(__file__) is True

    def test_exist_ok_not_strict(self):
        """
        测试非严格模式下，文件路径存在
        """
        use_not_strict_mode()
        assert assert_exist(__file__) is True

    def test_exist_error_strict(self):
        """
        测试严格模式下，文件路径不存在
        """
        use_strict_mode()
        with pytest.raises(FileNotFoundError):
            assert_exist(generate())
        assert assert_exist(generate(), ra=False) is False

    def test_exist_error_not_strict(self):
        """
        测试非严格模式下，文件路径不存在
        """
        use_not_strict_mode()
        assert assert_exist(generate()) is False

    def test_file_error_strict(self):
        """
        测试严格模式下，文件路径是文件夹
        """
        use_strict_mode()
        with pytest.raises(NotADirectoryError):
            assert_exist(__file__, target_type="folder")

    def test_file_error_not_strict(self):
        """
        测试非严格模式下，文件路径是文件夹
        """
        use_not_strict_mode()
        assert assert_exist(__file__, target_type="folder") is True

    def test_folder_error_strict(self):
        """
        测试严格模式下，文件夹路径是文件
        """
        use_strict_mode()
        with pytest.raises(IsADirectoryError):
            assert_exist(os.path.dirname(__file__), target_type="file")

    def test_folder_error_not_strict(self):
        """
        测试非严格模式下，文件夹路径是文件
        """
        use_not_strict_mode()
        assert assert_exist(os.path.dirname(__file__), target_type="file") is True
