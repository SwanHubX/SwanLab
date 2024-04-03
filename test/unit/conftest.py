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
import os
from tutils.config import TEMP_PATH
import shutil


@pytest.fixture(scope="session", autouse=True)
def setup_before_all():
    # 在整个 pytest 运行之前执行的前置操作
    if os.path.exists(TEMP_PATH):
        shutil.rmtree(TEMP_PATH)
    os.mkdir(TEMP_PATH)
