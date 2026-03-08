"""
@author: cunyue
@file: test_init.py
@time: 2026/3/5
@description: 测试 swanlab.sdk.__init__ 模块中的工具函数
"""

from swanlab.sdk.cmd.init import set_nested_value


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
