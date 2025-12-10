#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:44
@File: pytest_text.py
@IDE: pycharm
@Description:
    测试文本处理模块
"""
import pytest

from swanlab.data.modules import Text


def test_text_ok():
    # ---------------------------------- 字符串输入 ----------------------------------
    mock = "这是一段测试文本"
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == mock
    assert buffer is None
    assert text.get_more() is None
    # ---------------------------------- float输入 ----------------------------------
    mock = 1.0
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == "1.0"
    assert buffer is None
    assert text.get_more() is None
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock)
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more() is None


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
    # ---------------------------------- float输入 ----------------------------------
    mock = 1.0
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1.0"
    assert buffer is None
    assert text.get_more()["caption"] == "test"
    # ---------------------------------- int输入 ----------------------------------
    mock = 1
    text = Text(data=mock, caption="test")
    data, buffer = text.parse()
    assert data == "1"
    assert buffer is None
    assert text.get_more()["caption"] == "test"


def test_text_nested():
    """测试Text类支持嵌套输入（套娃）"""
    # 创建基础Text实例
    base_text = Text(data="Hello World", caption="original")
    
    # 测试嵌套输入
    nested_text = Text(base_text)
    
    # 验证属性被正确复制
    assert nested_text.text_data == base_text.text_data
    assert nested_text.caption == base_text.caption
    
    # 测试可以覆盖caption
    nested_with_new_caption = Text(base_text, caption="new caption")
    assert nested_with_new_caption.text_data == base_text.text_data
    assert nested_with_new_caption.caption == "new caption"
    
    # 验证原始实例未被修改
    assert base_text.caption == "original"
    
    # 测试数字类型的嵌套
    base_num = Text(data=42, caption="number")
    nested_num = Text(base_num)
    assert nested_num.text_data == "42"
    assert nested_num.caption == "number"
