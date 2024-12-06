#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/9/6 16:34
@File: test_cli_task.py
@IDE: pycharm
@Description:
    测试 swanlab task 相关函数
"""
import pytest

import tutils
from swanlab.cli.commands.task.utils import OutputModel


@pytest.mark.skipif(tutils.is_skip_cloud_test, reason="skip cloud test")
def test_output_model_none():
    om = OutputModel("123456", {})
    assert om.cuid == "123456"
    assert om.path is None
    assert om.size is None
    assert om.output_url is None


@pytest.mark.skipif(tutils.is_skip_cloud_test, reason="skip cloud test")
def test_output_model_ok():
    om = OutputModel("123456", {"path": "nothing.zip", "size": 123})
    assert om.cuid == "123456"
    assert om.path == "nothing.zip"
    assert om.size == OutputModel.fmt_size(123)


@pytest.mark.skipif(tutils.is_skip_cloud_test, reason="skip cloud test")
def test_output_model_fmt_size():
    assert OutputModel.fmt_size(1) == "1.00 Byte"
    assert OutputModel.fmt_size(1024) == "1.00 KB"
    assert OutputModel.fmt_size(1024 * 1024) == "1.00 MB"
    assert OutputModel.fmt_size(1024 * 1024 * 1024) == "1.00 GB"
    assert OutputModel.fmt_size(1024 * 1024 * 1024 * 1024) == "1.00 TB"
