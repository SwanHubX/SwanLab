"""
@author: cunyue
@file: test_scope.py
@time: 2026/3/8
@description: 测试 Scope 上下文作用域管理器的隔离、嵌套与冒泡机制
"""

from swanlab.sdk.internal.pkg.scope import (
    Scope,
    get_context,
    set_context,
)


def test_no_scope_safe_fallback():
    """测试在没有任何 Scope 的情况下调用，是否会安全地静默忽略/返回默认值"""
    # 写入不该报错
    set_context("ghost_key", "ghost_value")

    # 读取应该返回 None 或指定的默认值
    assert get_context("ghost_key") is None
    assert get_context("ghost_key", default="safe") == "safe"


def test_basic_scope():
    """测试最基础的作用域读写"""
    with Scope() as scope:
        set_context("trace_id", "12345")
        set_context("user", "cunyue")

        # 验证内部直接读取
        assert get_context("trace_id") == "12345"
        # 验证 scope 实例的数据字典
        assert scope.data == {"trace_id": "12345", "user": "cunyue"}

    # 退出 with 块后，上下文应该被清空
    assert get_context("trace_id") is None


def test_nested_scope_isolation_and_lookup():
    """测试嵌套作用域的隔离性（默认不冒泡）以及向上追溯读取（Chain Lookup）"""
    with Scope() as outer:
        set_context("env", "production")
        set_context("step", "outer_step")

        with Scope() as inner:
            set_context("step", "inner_step")
            set_context("inner_only", "secret")

            # 1. 内部同名 key 会覆盖外部的 key（作用域遮蔽）
            assert get_context("step") == "inner_step"

            # 2. 内部没有的 key，会自动向上追溯找到 outer 的 key
            assert get_context("env") == "production"

            # 3. 验证 inner 的实际 data 字典只存了自己的一份
            assert inner.data == {"step": "inner_step", "inner_only": "secret"}

        # 4. 退出 inner 后，outer 的数据保持原样，没有被 inner 污染
        assert get_context("step") == "outer_step"
        assert get_context("inner_only") is None
        assert outer.data == {"env": "production", "step": "outer_step"}


def test_nested_scope_bubble_up():
    """测试嵌套作用域开启 bubble_up=True 时的冒泡合并机制"""
    with Scope() as outer:
        set_context("base_info", "swanlab")

        with Scope(bubble_up=True):
            set_context("deep_info", "metrics")
            set_context("base_info", "overridden_by_inner")  # 测试键冲突覆盖

        # 退出 inner 后，由于开启了 bubble_up，inner 的数据应该合并到了 outer
        assert outer.data["deep_info"] == "metrics"
        # 冲突的 key 会以后退出的 inner 为准
        assert outer.data["base_info"] == "overridden_by_inner"


def test_deep_function_call_simulation():
    """模拟真实的深层函数调用，验证‘透传污染’问题的解决效果"""

    def deep_level_3():
        # 在极深处隔空上报数据
        set_context("db_query_time", 0.45)
        set_context("warning", "slow query")
        return "data_from_db"

    def mid_level_2():
        return deep_level_3()

    def top_level_1():
        return mid_level_2()

    # 外层业务代码
    with Scope() as scope:
        result = top_level_1()

        # 验证正常的返回值不受影响
        assert result == "data_from_db"

        # 验证隔空取物成功
        assert get_context("db_query_time") == 0.45
        assert get_context("warning") == "slow query"
        assert scope.data == {"db_query_time": 0.45, "warning": "slow query"}
