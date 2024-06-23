#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 16:52
@File: conftest.py.py
@IDE: pycharm
@Description:
    配置pytest
"""
import pytest
from tutils import TEMP_PATH, reset_some_env
import swanlab.env as E
import shutil
import os


def count_files_in_directory(directory, exclude_prefixes=("__", ".", "~")):
    """
    计算目录下的文件数量
    """
    file_count, folder_count = 0, 0
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if any(filename.startswith(prefix) for prefix in exclude_prefixes):
            continue
        if os.path.isfile(file_path):
            file_count += 1
        elif os.path.isdir(file_path):
            folder_count += 1
            r = count_files_in_directory(file_path)
            file_count += r[0]
            folder_count += r[1]
    return file_count, folder_count


@pytest.fixture(scope="session", autouse=True)
def setup():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # 记住当前文件夹下的所有文件数量，用于测试结束后比较，确保测试结束后没有多余文件
    pre_file_count, pre_folder_count = count_files_in_directory(current_dir)
    yield
    # 确保测试结束后没有多余文件
    now_file_count, now_folder_count = count_files_in_directory(current_dir)
    assert pre_file_count == now_file_count
    assert pre_folder_count == now_folder_count


@pytest.fixture(scope="function", autouse=True)
def setup_each():
    """
    对每一个测试函数进行设置
    对于每一个测试函数，将环境变量恢复原状，清空temp文件夹后重新创建一个新的
    """
    for key in E.SwanLabEnv.list():
        if key in os.environ:
            del os.environ[key]
    reset_some_env()
    # 清空temp文件夹
    if os.path.exists(TEMP_PATH):
        shutil.rmtree(TEMP_PATH)
    os.mkdir(TEMP_PATH)
    yield
