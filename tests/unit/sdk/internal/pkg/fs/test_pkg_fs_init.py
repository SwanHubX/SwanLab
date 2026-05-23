"""
@author: cunyue
@file: test_pkg_fs_init.py
@time: 2026/5/23
@description: 测试 fs.__init__ 中的 safe_fmt 和 safe_truncate
"""

import pytest

from swanlab.sdk.internal.pkg.fs import safe_fmt, safe_truncate


class TestSafeFmt:
    def test_ascii_passthrough(self):
        assert safe_fmt("hello") == "hello"

    def test_illegal_chars_replaced(self):
        result = safe_fmt("a<b>c:d")
        assert result == "a_b_c_d"

    def test_dot_and_dash_replaced(self):
        result = safe_fmt("my-dir.name")
        assert result == "my_dir_name"

    def test_consecutive_underscores_collapsed(self):
        result = safe_fmt("a<>b<>c")
        assert result == "a_b_c"

    def test_leading_trailing_stripped(self):
        result = safe_fmt("---hello---")
        assert result == "hello"

    def test_chinese_preserved(self):
        result = safe_fmt("我的电脑")
        assert result == "我的电脑"

    def test_all_illegal_uses_fallback(self):
        result = safe_fmt("<>:|?*#%", fallback="my_fb")
        assert result == "my_fb"

    def test_all_illegal_default_fallback(self):
        result = safe_fmt("<>:|?*#%")
        assert result == "unknown"

    def test_custom_fallback(self):
        result = safe_fmt("<>:|?*#%", fallback="unknown_host")
        assert result == "unknown_host"

    def test_empty_string_uses_fallback(self):
        assert safe_fmt("") == "unknown"

    def test_whitespace_only_uses_fallback(self):
        assert safe_fmt("   ") == "unknown"

    def test_mixed_chinese_and_illegal(self):
        result = safe_fmt("我的:电脑/主机")
        assert result == "我的_电脑_主机"


class TestSafeTruncate:
    def test_short_name_unchanged(self):
        result, truncated = safe_truncate("hello", 255)
        assert result == "hello"
        assert truncated is False

    def test_exact_length_unchanged(self):
        name = "a" * 10
        result, truncated = safe_truncate(name, 10)
        assert result == name
        assert truncated is False

    def test_truncate_basic(self):
        name = "a" * 20
        result, truncated = safe_truncate(name, 10)
        assert truncated is True
        assert result.startswith("aaa...")
        assert len(result.encode("utf-8")) <= 10

    def test_truncate_format(self):
        name = "abcdefghijklmnopqrstuvwxyz"
        result, truncated = safe_truncate(name, 10)
        assert truncated is True
        assert "..." in result
        assert result.startswith("abc")

    def test_max_length_too_small(self):
        with pytest.raises(ValueError, match="must be >= 7"):
            safe_truncate("a" * 20, 5)

    def test_max_length_exactly_7(self):
        name = "a" * 20
        result, truncated = safe_truncate(name, 7)
        assert truncated is True
        assert len(result.encode("utf-8")) <= 7
        assert "..." in result

    def test_chinese_within_limit(self):
        name = "测试" * 100
        result, truncated = safe_truncate(name, 30)
        assert truncated is True
        assert len(result.encode("utf-8")) <= 30

    def test_chinese_prefix_exceeds_max_length(self):
        name = "测试" * 100
        with pytest.raises(ValueError, match="first 3 characters"):
            safe_truncate(name, 7)

    def test_byte_length_not_character_length(self):
        name = "测试" * 50
        result, truncated = safe_truncate(name, 20)
        assert truncated is True
        assert len(result.encode("utf-8")) <= 20
        assert len(result) != len(result.encode("utf-8"))

    def test_default_max_length(self):
        name = "a" * 300
        result, truncated = safe_truncate(name)
        assert truncated is True
        assert len(result.encode("utf-8")) <= 255
