"""
@author: cunyue
@file: test_run_helper.py
@time: 2026/3/12 16:54
@description: 测试run工具函数
"""

from unittest.mock import patch

import pytest

import swanlab.sdk.internal.run.utils_fmt as fmt

# ─────────────────────────── flatten_dict ────────────────────────────


class TestFlattenDict:
    """测试 fmt.flatten_dict 函数"""

    def test_empty_dict(self):
        """空字典直接返回空字典"""
        assert fmt.flatten_dict({}) == {}

    def test_flat_dict_unchanged(self):
        """无嵌套字典不做任何变换"""
        assert fmt.flatten_dict({"a": 1, "b": "x"}) == {"a": 1, "b": "x"}

    def test_nested_one_level(self):
        """单层嵌套展开"""
        assert fmt.flatten_dict({"a": {"b": 1}}) == {"a/b": 1}

    def test_nested_deep(self):
        """深层嵌套展开"""
        assert fmt.flatten_dict({"a": {"b": {"c": 1}}}) == {"a/b/c": 1}

    def test_mixed_flat_and_nested(self):
        """平铺与嵌套混合"""
        result = fmt.flatten_dict({"x": 0, "a": {"b": 1, "c": 2}})
        assert result == {"x": 0, "a/b": 1, "a/c": 2}

    def test_int_key_converted_to_str(self):
        """整数 key 自动转为字符串"""
        assert fmt.flatten_dict({1: {2: "v"}}) == {"1/2": "v"}  # type: ignore

    def test_duplicate_key_latter_wins(self):
        """出现同名展开键时，后者覆盖前者，并触发一次警告"""
        # 先展开嵌套得到 a/b=1，再直接赋值 "a/b"=99 → 发生覆盖
        data = {"a": {"b": 1}, "a/b": 99}
        with patch("swanlab.sdk.internal.run.utils_fmt.console") as mock_console:
            result = fmt.flatten_dict(data)
        mock_console.warning.assert_called_once()
        assert result["a/b"] == 99

    def test_returns_parent_dict_when_provided(self):
        """传入 parent_dict 时结果合并到其中"""
        parent = {"existing": "value"}
        result = fmt.flatten_dict({"a": 1}, parent_dict=parent)
        assert result is parent
        assert result == {"existing": "value", "a": 1}


# ─────────────────────────── validate_key ────────────────────────────


class TestValidateKey:
    """测试 fmt.validate_key 函数"""

    @pytest.fixture(autouse=True)
    def reset_warned_keys(self):
        """每个测试前清空全局警告缓存，避免用例间相互污染"""
        fmt._WARNED_KEYS.clear()
        yield
        fmt._WARNED_KEYS.clear()

    def test_valid_key_unchanged(self):
        """合法 key 原样返回"""
        assert fmt.validate_key("loss") == "loss"
        assert fmt.validate_key("train/acc") == "train/acc"
        assert fmt.validate_key("v1.0_metric") == "v1.0_metric"

    def test_strip_leading_trailing(self):
        """头尾的空白、. 和 / 被剥离"""
        with patch("swanlab.sdk.internal.run.utils_fmt.console"):
            assert fmt.validate_key("  loss  ") == "loss"
            assert fmt.validate_key("./key/") == "key"
            assert fmt.validate_key("...key...") == "key"

    @pytest.mark.parametrize(
        "invalid, expected",
        [
            ("my key", "my_key"),  # 空格 -> _
            ("a@b#c", "a_b_c"),  # 特殊符号 -> _
        ],
    )
    def test_invalid_chars_replaced(self, invalid, expected):
        """非法字符被替换为下划线"""
        with patch("swanlab.sdk.internal.run.utils_fmt.console"):
            result = fmt.validate_key(invalid)
        assert result == expected

    def test_unicode_word_chars_allowed(self):
        """Python \\w 匹配 Unicode 字符（汉字等），不会被替换"""
        # \w 在 Python3 默认模式下匹配 Unicode 字母数字，汉字合法
        result = fmt.validate_key("中文key")
        assert result == "中文key"

    def test_truncate_long_key(self):
        """超长 key 被截断到 max_len"""
        long_key = "a" * 300
        result = fmt.validate_key(long_key, max_len=255)
        assert len(result) == 255
        assert result == "a" * 255

    def test_custom_max_len(self):
        """自定义 max_len 生效"""
        with patch("swanlab.sdk.internal.run.utils_fmt.console"):
            result = fmt.validate_key("abcdef", max_len=3)
        assert result == "abc"

    def test_non_string_input_int(self):
        """整数 key 被转为字符串处理"""
        assert fmt.validate_key(42) == "42"  # type: ignore

    def test_non_string_input_float(self):
        """浮点 key 被转为字符串处理"""
        assert fmt.validate_key(3.14) == "3.14"  # type: ignore

    def test_empty_after_sanitization_raises(self):
        """清洗后为空时抛出 ValueError"""
        with pytest.raises(ValueError, match="invalid or empty after sanitization"):
            fmt.validate_key("   ")
        with pytest.raises(ValueError, match="invalid or empty after sanitization"):
            fmt.validate_key("...")

    def test_warning_fires_once_per_key(self):
        """相同非法 key 只触发一次警告"""
        with patch("swanlab.sdk.internal.run.utils_fmt.console") as mock_console:
            fmt.validate_key("bad key")
            fmt.validate_key("bad key")
            fmt.validate_key("bad key")
        assert mock_console.warning.call_count == 1

    def test_warning_fires_for_different_keys(self):
        """不同非法 key 各自独立触发一次警告"""
        with patch("swanlab.sdk.internal.run.utils_fmt.console") as mock_console:
            fmt.validate_key("bad key1")
            fmt.validate_key("bad key2")
        assert mock_console.warning.call_count == 2
