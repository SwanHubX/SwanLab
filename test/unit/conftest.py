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
from tutils import clear, init_db, SWANLAB_DIR, SWANLAB_LOG_DIR, PACKAGE_PATH, TEMP_PATH
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
        if os.path.isfile(file_path):
            if any(filename.startswith(prefix) for prefix in exclude_prefixes):
                continue
            file_count += 1
        elif os.path.isdir(file_path):
            folder_count += 1
            r = count_files_in_directory(file_path)
            file_count += r[0]
            folder_count += r[1]
    return file_count, folder_count


@pytest.fixture(scope="session", autouse=True)
def setup_before_all():
    clear()
    init_db()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # 记住当前文件夹下的所有文件数量，用于测试结束后比较，确保测试结束后没有多余文件
    pre_file_count, pre_folder_count = count_files_in_directory(current_dir)
    yield
    # 确保测试结束后没有多余文件
    now_file_count, now_folder_count = count_files_in_directory(current_dir)
    assert pre_file_count == now_file_count
    assert pre_folder_count == now_folder_count


@pytest.fixture(scope="function", autouse=True)
def setup_before_each():
    E.reset_env()
    if os.path.exists(SWANLAB_DIR):
        shutil.rmtree(SWANLAB_DIR)
    yield
    E.reset_env()
    shutil.rmtree(SWANLAB_DIR, ignore_errors=True)
    os.environ[E.DEV] = "TRUE"
    os.environ[E.ROOT] = SWANLAB_LOG_DIR
    os.environ[E.PACKAGE] = PACKAGE_PATH
    os.environ[E.HOME] = TEMP_PATH
