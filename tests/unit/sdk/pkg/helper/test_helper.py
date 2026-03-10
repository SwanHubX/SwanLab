"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 20:47
@description: SwanLab SDK 辅助函数
"""

from swanlab.sdk.utils.helper import strip_none


class TestStripNone:
    """测试 strip_none 函数"""

    def test_empty_dict(self):
        """测试空字典"""
        assert strip_none({}) == {}

    def test_simple_none_values(self):
        """测试简单的 None 值会被移除"""
        data = {"a": 1, "b": None, "c": "test"}
        assert strip_none(data) == {"a": 1, "c": "test"}

    def test_nested_dict_with_none(self):
        """测试嵌套字典中的 None 值和空字典会被移除"""
        data = {"a": 1, "b": {"x": None, "y": 2}, "c": {"z": None}, "d": None}
        assert strip_none(data) == {"a": 1, "b": {"y": 2}}

    def test_deeply_nested_dict(self):
        """测试深层嵌套的字典"""
        data = {"level1": {"level2": {"level3": {"a": None, "b": 1}, "c": None}, "d": 2}}
        expected = {"level1": {"level2": {"level3": {"b": 1}}, "d": 2}}
        assert strip_none(data) == expected

    def test_preserve_falsy_values(self):
        """测试保留 False、0、空字符串等假值（它们不是 None）"""
        data = {"false": False, "zero": 0, "empty_str": "", "none": None}
        assert strip_none(data) == {"false": False, "zero": 0, "empty_str": ""}
