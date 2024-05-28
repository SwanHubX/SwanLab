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
from tutils import clear, init_db, SWANLAB_DIR
from swanlab.env import reset_env
import shutil
import os


@pytest.fixture(scope="session", autouse=True)
def setup_before_all():
    clear()
    init_db()


@pytest.fixture(scope="function", autouse=True)
def setup_before_each():
    reset_env()
    if os.path.exists(SWANLAB_DIR):
        shutil.rmtree(SWANLAB_DIR)
    yield
    reset_env()
    shutil.rmtree(SWANLAB_DIR, ignore_errors=True)
