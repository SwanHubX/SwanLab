#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/28 16:21
@File: pytest_file.py
@IDE: pycharm
@Description:
    测试文件格式函数
"""
import pytest
from nanoid import generate
from swanlab.utils.file import check_proj_name_format, check_exp_name_format, _auto_cut


class TestAutoCut:
    @pytest.mark.parametrize("name, value", [
        [generate(), generate(size=101)],
        [generate(), generate(size=1000)],
        [generate(), generate(size=10000)],
    ])
    def test_cut(self, name: str, value: str):
        """
        测试自动截断
        """
        assert len(_auto_cut(name, value, 100, True)) == 100

    @pytest.mark.parametrize("name, value", [
        [generate(), generate(size=101)],
        [generate(), generate(size=1000)],
        [generate(), generate(size=10000)],
    ])
    def test_no_cut(self, name: str, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError) as e:
            _auto_cut(name, value, 100, False)
        assert name in str(e.value)


class TestProjName:
    @pytest.mark.parametrize("value", [
        generate(size=100),
        generate(size=1),
        "-",
        "_",
    ])
    def test_proj_name_common(self, value):
        """
        测试正常情况
        """
        assert check_proj_name_format(value) == value

    @pytest.mark.parametrize("value", [
        None,
        1,
        [],
        {}
    ])
    def test_proj_name_type_error(self, value: str):
        """
        测试类型错误
        """
        with pytest.raises(TypeError):
            check_proj_name_format(value)

    @pytest.mark.parametrize("value", [
        "",
        "啊哈哈",
        "&^%",
        "/;]x]"
    ])
    def test_proj_name_value_error(self, value: str):
        """
        测试空值或者不合法值
        """
        with pytest.raises(ValueError):
            check_proj_name_format(value)

    @pytest.mark.parametrize("value", [
        generate(size=101),
        generate(size=1000),
        generate(size=10000),
    ])
    def test_proj_name_auto_cut(self, value: str):
        """
        测试自动截断
        """
        assert len(check_proj_name_format(value)) == 100

    @pytest.mark.parametrize("value", [
        generate(size=101),
        generate(size=1000),
        generate(size=10000),
    ])
    def test_proj_name_no_cut(self, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError):
            check_proj_name_format(value, auto_cut=False)
