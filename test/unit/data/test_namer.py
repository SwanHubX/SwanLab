#!/usr/bin/env python
r"""
@DATE: 2024/9/20 14:42
@File: test_namer.py
@IDE: pycharm
@Description:
    测试命名器、取色器
"""

import pytest

from swanlab.data import namer
from swanlab.data.namer import hex_to_rgb


def test_name_no_index():
    name = namer.generate_name()
    assert isinstance(name, str)


def test_name_with_index():
    name = namer.generate_name(4)
    assert isinstance(name, str)
    # 极大数
    name = namer.generate_name(999999999)
    assert isinstance(name, str)


def test_color_no_index():
    colors = namer.generate_colors()
    assert len(colors) == 2
    assert isinstance(colors, tuple)


def test_color_with_index():
    colors = namer.generate_colors(4)
    assert len(colors) == 2
    assert isinstance(colors, tuple)
    # 极大数
    colors = namer.generate_colors(999999999)
    assert len(colors) == 2
    assert isinstance(colors, tuple)


def test_hex_to_rgb_normal():
    """测试正常的hex颜色转换"""
    # 测试标准6位十六进制
    assert hex_to_rgb("#528d59") == (82, 141, 89)
    assert hex_to_rgb("528d59") == (82, 141, 89)

    # 测试3位简写形式
    assert hex_to_rgb("#fff") == (255, 255, 255)
    assert hex_to_rgb("fff") == (255, 255, 255)

    # 测试黑白色
    assert hex_to_rgb("#000000") == (0, 0, 0)
    assert hex_to_rgb("#FFFFFF") == (255, 255, 255)

    # 测试大小写混合
    assert hex_to_rgb("#fF0000") == (255, 0, 0)


def test_hex_to_rgb_edge_cases():
    """测试边界情况"""
    # 测试带空格的情况
    assert hex_to_rgb(" #528d59 ") == (82, 141, 89)
    assert hex_to_rgb(" 528d59 ") == (82, 141, 89)


def test_hex_to_rgb_invalid_input():
    """测试无效输入"""
    # 测试长度错误
    with pytest.raises(ValueError, match="Invalid hex color length"):
        hex_to_rgb("#12345")  # 5位

    with pytest.raises(ValueError, match="Invalid hex color length"):
        hex_to_rgb("1234567")  # 7位

    # 测试非法字符
    with pytest.raises(ValueError, match="Invalid hex color"):
        hex_to_rgb("#xyz123")

    with pytest.raises(ValueError, match="Invalid hex color"):
        hex_to_rgb("#12345g")

    # 测试空字符串
    with pytest.raises(ValueError, match="Invalid hex color length"):
        hex_to_rgb("")

    # 测试None
    with pytest.raises(AttributeError):
        hex_to_rgb(None)


def test_hex_to_rgb_rgb_bounds():
    """测试RGB值的范围"""
    # 测试最小值(0)
    result = hex_to_rgb("#000000")
    assert all(0 <= v <= 255 for v in result)

    # 测试最大值(255)
    result = hex_to_rgb("#FFFFFF")
    assert all(0 <= v <= 255 for v in result)

    # 测试中间值
    result = hex_to_rgb("#7F7F7F")
    assert all(0 <= v <= 255 for v in result)


def test_hex_to_rgb_basic_colors():
    """测试基本颜色"""
    assert hex_to_rgb("#FF0000") == (255, 0, 0)  # 红
    assert hex_to_rgb("#00FF00") == (0, 255, 0)  # 绿
    assert hex_to_rgb("#0000FF") == (0, 0, 255)  # 蓝
    assert hex_to_rgb("#FFFF00") == (255, 255, 0)  # 黄
    assert hex_to_rgb("#FF00FF") == (255, 0, 255)  # 品红
    assert hex_to_rgb("#00FFFF") == (0, 255, 255)  # 青


class TestGenerateRunId:

    def test_run_id_length_is_21_characters(self):
        run_id = namer.generate_run_id()
        assert len(run_id) == 21

    def test_run_id_contains_only_lowercase_and_digits(self):
        run_id = namer.generate_run_id()
        assert all(c.islower() or c.isdigit() for c in run_id)

    def test_run_id_is_random_each_time(self):
        run_id_1 = namer.generate_run_id()
        run_id_2 = namer.generate_run_id()
        assert run_id_1 != run_id_2
