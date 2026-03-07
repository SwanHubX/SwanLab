"""
@author: cunyue
@file: test_init.py
@time: 2026/3/5
@description: 测试 swanlab.sdk.__init__ 模块中的工具函数
"""

from swanlab.sdk.cmd.init import set_nested_value, strip_none


class TestSetNestedValue:
    """测试 set_nested_value 函数"""

    def test_simple_key(self):
        """测试简单的单层键"""
        d = {}
        set_nested_value(d, "a", 1)
        assert d == {"a": 1}

    def test_nested_key(self):
        """测试嵌套键路径"""
        d = {}
        set_nested_value(d, "a.b.c", 123)
        assert d == {"a": {"b": {"c": 123}}}

    def test_update_existing_nested_value(self):
        """测试更新已存在的嵌套值"""
        d = {"a": {"b": {"c": 1}}}
        set_nested_value(d, "a.b.c", 999)
        assert d == {"a": {"b": {"c": 999}}}

    def test_partial_existing_path(self):
        """测试部分路径已存在的情况"""
        d = {"a": {"x": 1}}
        set_nested_value(d, "a.b.c", 2)
        assert d == {"a": {"x": 1, "b": {"c": 2}}}

    def test_none_value_not_set(self):
        """测试 None 值不会被设置"""
        d = {}
        set_nested_value(d, "a", None)
        assert d == {}

    def test_none_value_nested_not_set(self):
        """测试嵌套路径的 None 值不会被设置"""
        d = {"x": 1}
        set_nested_value(d, "a.b.c", None)
        assert d == {"x": 1}


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
