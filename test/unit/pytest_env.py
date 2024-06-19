#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/24 20:25
@File: pytest_env.py
@IDE: pycharm
@Description:
    测试swanlab/env.py中的环境变量和工具函数
"""
import pytest
import swanlab.env as E
from nanoid import generate
import tutils as T
import os


class TestPackagePath:
    def test_default(self):
        """
        测试默认的package路径
        """
        del os.environ[E.SwanLabEnv.SWANLAB_PACKAGE.value]
        project_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        package_path = os.path.join(project_dir, "swanlab", "package.json")
        assert E.get_package_path() == package_path

    def test_use_env(self):
        """
        使用环境变量指定package路径
        """
        os.environ[E.SwanLabEnv.SWANLAB_PACKAGE.value] = T.PACKAGE_PATH

    def test_path_is_dir(self):
        """
        package路径是一个目录
        """
        os.environ[E.SwanLabEnv.SWANLAB_PACKAGE.value] = T.TEMP_PATH
        with pytest.raises(IsADirectoryError):
            E.get_package_path()

    def test_path_not_exists(self):
        """
        package路径不存在
        """
        os.environ[E.SwanLabEnv.SWANLAB_PACKAGE.value] = generate()
        with pytest.raises(FileNotFoundError):
            E.get_package_path()
