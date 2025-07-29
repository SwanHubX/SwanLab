#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/19 14:43
@File: pytest_formater.py
@IDE: pycharm
@Description:

"""
import pytest
from nanoid import generate

# noinspection PyProtectedMember
from swanlab.formatter import (
    check_proj_name_format,
    _auto_cut,
    check_key_format,
    check_exp_name_format,
    check_run_id_format,
)


class TestAutoCut:
    @pytest.mark.parametrize(
        "name, value",
        [
            [generate(), generate(size=101)],
            [generate(), generate(size=1000)],
            [generate(), generate(size=10000)],
        ],
    )
    def test_cut(self, name: str, value: str):
        """
        测试自动截断
        """
        assert len(_auto_cut(name, value, 100, True)) == 100

    @pytest.mark.parametrize(
        "name, value",
        [
            [generate(), generate(size=101)],
            [generate(), generate(size=1000)],
            [generate(), generate(size=10000)],
        ],
    )
    def test_no_cut(self, name: str, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError) as e:
            _auto_cut(name, value, 100, False)
        assert name in str(e.value)


class TestProjName:
    @pytest.mark.parametrize(
        "value",
        [generate(size=100), generate(size=1), "-", "_", ".12", "1", "1.b", "a.b", "+", "1+1"],
    )
    def test_proj_name_common(self, value):
        """
        测试正常情况
        """
        assert check_proj_name_format(value) == value

    @pytest.mark.parametrize("value", [None, 1, [], {}])
    def test_proj_name_type_error(self, value: str):
        """
        测试类型错误
        """
        with pytest.raises(TypeError):
            check_proj_name_format(value)

    @pytest.mark.parametrize("value", ["", " " "啊哈哈", "&^%", "/;]x]"])
    def test_proj_name_value_error(self, value: str):
        """
        测试空值或者不合法值
        """
        with pytest.raises(ValueError):
            check_proj_name_format(value)

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=101),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_proj_name_auto_cut(self, value: str):
        """
        测试自动截断
        """
        assert len(check_proj_name_format(value)) == 100

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=101),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_proj_name_no_cut(self, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError):
            check_proj_name_format(value, auto_cut=False)


class TestExpName:
    @pytest.mark.parametrize(
        "value",
        [
            generate(size=250),
            generate(size=95),
            generate(size=1),
            "-",
            "_",
            ".12",
            "1",
            "1.b",
            "a.b",
            "+",
            "1+1",
            "你好",
        ],
    )
    def test_exp_name_common(self, value):
        """
        测试正常情况
        """
        assert check_exp_name_format(value) == value

    @pytest.mark.parametrize("value", [None, 1, [], {}])
    def test_exp_name_type_error(self, value: str):
        """
        测试类型错误
        """
        with pytest.raises(TypeError):
            check_exp_name_format(value)

    @pytest.mark.parametrize("value", ["", " "])
    def test_exp_name_value_error(self, value: str):
        """
        测试空值
        """
        with pytest.raises(ValueError):
            check_exp_name_format(value)

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=251),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_exp_name_auto_cut(self, value: str):
        """
        测试自动截断
        """
        assert len(check_exp_name_format(value)) == 250

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=251),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_exp_name_no_cut(self, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError):
            check_exp_name_format(value, auto_cut=False)


class TestTag:

    @pytest.mark.parametrize(
        "value", [generate(size=255), generate(size=100), generate(size=1), "12", "-", "_", "👾👾👾👾👾👾"]
    )
    def test_tag_common(self, value):
        """
        测试正常情况
        """
        assert check_key_format(value) == value

    def test_key_blank(self):
        """
        测试收尾空格情况
        """
        with pytest.raises(ValueError):
            check_key_format("  ")
        key = "  " + "abc" + "  "
        assert check_key_format(key) == "abc"

    @pytest.mark.parametrize(
        "value",
        [
            None,
            1,
            [],
            {},
        ],
    )
    def test_tag_type_error(self, value: str):
        """
        测试类型错误
        """
        with pytest.raises(TypeError):
            check_key_format(value)

    @pytest.mark.parametrize("value", ["", "   ", " " * 256, ".sas", "/asa", "abc/", "bac."])
    def test_tag_value_error(self, value: str):
        """
        测试不合法值
        """
        with pytest.raises(ValueError):
            check_key_format(value)

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=256),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_tag_auto_cut(self, value: str):
        """
        测试自动截断
        """
        assert len(check_key_format(value)) == 255

    @pytest.mark.parametrize(
        "value",
        [
            generate(size=256),
            generate(size=1000),
            generate(size=10000),
        ],
    )
    def test_tag_no_cut(self, value: str):
        """
        测试不自动截断
        """
        with pytest.raises(IndexError):
            check_key_format(value, auto_cut=False)


class TestRunIdFormat:
    @staticmethod
    def test_run_id_format_valid_run_id():
        assert check_run_id_format("abcdefghijklmnopqrstu") == "abcdefghijklmnopqrstu"

    @staticmethod
    def test_run_id_format_invalid_length():
        with pytest.raises(
            ValueError, match=r"id .* is invalid, it must be 21 characters of lowercase letters and digits"
        ):
            check_run_id_format("shortid")

    @staticmethod
    def test_run_id_format_invalid_characters():
        with pytest.raises(
            ValueError, match=r"id .* is invalid, it must be 21 characters of lowercase letters and digits"
        ):
            check_run_id_format("abc123!@#def456ghi789")

    @staticmethod
    def test_run_id_format_none_input():
        assert check_run_id_format(None) is None

    @staticmethod
    def test_run_id_format_numeric_input():
        assert check_run_id_format(123456789012345678901) == "123456789012345678901"

    @staticmethod
    def test_run_id_format_empty_string():
        with pytest.raises(
            ValueError, match=r"id .* is invalid, it must be 21 characters of lowercase letters and digits"
        ):
            check_run_id_format("")
