#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:44
@File: pytest_text.py
@IDE: pycharm
@Description:
    测试文本处理模块
"""
from swanlab.data.modules import Text
import pytest


def test_text_ok():
    # ---------------------------------- 字符串输入 ----------------------------------
    mock = "这是一段测试文本"
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == mock
    assert buffer is None
    assert text.get_more() is None
    assert text.get_config() is None
    # ---------------------------------- float输入 ----------------------------------
    mock = 1.0
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == "1.0"
    assert buffer is None
    assert text.get_more() is None
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more() is None
    assert text.get_config() is None


def test_text_error_type():
    mock = [1, 2, 3]
    with pytest.raises(TypeError):
        Text(data=mock)  # noqa


def test_text_caption():
    mock = "这是一段测试文本"
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == mock
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- float输入 ----------------------------------
    mock = 1.0
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1.0"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    assert text.get_config() is None
