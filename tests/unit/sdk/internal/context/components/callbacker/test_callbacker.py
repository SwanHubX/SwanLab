"""
@author: cunyue
@file: test_callbacker.py
@time: 2026/4/30
@description: CallbackManager 及 create_callback_manager 单元测试
"""

import pytest

from swanlab.sdk.internal.context.components.callbacker import (
    CallbackManager,
    create_callback_manager,
    global_callbacker,
)
from swanlab.sdk.protocol import Callback


class _AlphaCallback(Callback):
    @property
    def name(self) -> str:
        return "alpha"


class _BetaCallback(Callback):
    @property
    def name(self) -> str:
        return "beta"


class _AlphaDupCallback(Callback):
    """与 Alpha 同名，用于测试覆盖"""

    @property
    def name(self) -> str:
        return "alpha"


class TestCallbackManager:
    """_CallbackManager 基础操作"""

    def test_merge_callbacks(self):
        mgr = CallbackManager()
        mgr.merge_callbacks([_AlphaCallback(), _BetaCallback()])
        assert {cb.name for cb in mgr.registered_callbacks} == {"alpha", "beta"}

    def test_merge_single_callback(self):
        """支持传入单个 Callback 对象"""
        mgr = CallbackManager()
        mgr.merge_callbacks(_AlphaCallback())
        assert len(mgr.registered_callbacks) == 1
        assert mgr.registered_callbacks[0].name == "alpha"

    def test_merge_overwrites_duplicate(self):
        mgr = CallbackManager()
        mgr.merge_callbacks([_AlphaCallback()])
        mgr.merge_callbacks([_AlphaDupCallback()])
        assert len(mgr.registered_callbacks) == 1
        assert isinstance(mgr.registered_callbacks[0], _AlphaDupCallback)

    def test_merge_none_or_empty(self):
        mgr = CallbackManager()
        mgr.merge_callbacks(None)
        mgr.merge_callbacks([])
        assert len(mgr.registered_callbacks) == 0

    def test_merge_invalid_type_raises(self):
        with pytest.raises(TypeError, match="Expected swanlab.Callback"):
            CallbackManager().merge_callbacks(["not_a_callback"])  # type: ignore

    def test_remove_callback(self):
        mgr = CallbackManager()
        mgr.merge_callbacks([_AlphaCallback()])
        mgr.remove_callback("alpha")
        assert len(mgr.registered_callbacks) == 0
        mgr.remove_callback("nonexistent")  # 不存在时不报错

    def test_dispatch_forwards_to_callback(self):
        calls = []

        class Rec(Callback):
            @property
            def name(self) -> str:
                return "rec"

            def on_run_finished(self, state, error=None):
                calls.append(state)

        mgr = CallbackManager()
        mgr.merge_callbacks([Rec()])
        mgr.on_run_finished("success")
        assert calls == ["success"]

    def test_dispatch_unknown_raises(self):
        with pytest.raises(AttributeError):
            CallbackManager().nonexistent_method()


class TestCreateCallbackManager:
    """create_callback_manager 工厂函数"""

    def test_empty_when_no_callbacks(self):
        assert len(create_callback_manager(None).registered_callbacks) == 0

    def test_local_only(self):
        mgr = create_callback_manager([_AlphaCallback()])
        assert {cb.name for cb in mgr.registered_callbacks} == {"alpha"}

    def test_local_single_callback(self):
        """支持传入单个 Callback 对象"""
        mgr = create_callback_manager(_AlphaCallback())
        assert {cb.name for cb in mgr.registered_callbacks} == {"alpha"}

    def test_global_and_local_merged(self):
        global_callbacker.merge_callbacks([_AlphaCallback()])
        try:
            mgr = create_callback_manager([_BetaCallback()])
            assert {cb.name for cb in mgr.registered_callbacks} == {"alpha", "beta"}
        finally:
            global_callbacker.remove_callback("alpha")

    def test_local_overwrites_global(self):
        global_callbacker.merge_callbacks([_AlphaCallback()])
        try:
            mgr = create_callback_manager([_AlphaDupCallback()])
            assert len(mgr.registered_callbacks) == 1
            assert isinstance(mgr.registered_callbacks[0], _AlphaDupCallback)
        finally:
            global_callbacker.remove_callback("alpha")

    def test_local_does_not_pollute_global(self):
        global_callbacker.merge_callbacks([_AlphaCallback()])
        try:
            create_callback_manager([_BetaCallback()])
            assert "beta" not in {cb.name for cb in global_callbacker.registered_callbacks}
        finally:
            global_callbacker.remove_callback("alpha")
