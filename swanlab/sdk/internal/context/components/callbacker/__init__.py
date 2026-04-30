"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:27
@description: SwanLab 回调模块
"""

from typing import TYPE_CHECKING, Dict, List, Optional

from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.protocol import Callback
from swanlab.sdk.typings.context import CallbacksType


class _CallbackManager:
    """
    SwanLab 回调函数的管理器。
    维护当前上下文中的全局回调状态，支持批量合并与指定移除。
    """

    def __init__(self):
        """
        初始化回调管理器
        """
        self._callbacks: Dict[str, Callback] = {}

    def merge_callbacks(self, callbacks: Optional[CallbacksType]) -> None:
        """批量合并回调函数到当前管理器中"""
        if callbacks is None:
            return
        if isinstance(callbacks, Callback):
            callbacks = [callbacks]

        for cb in callbacks:
            if not isinstance(cb, Callback):
                raise TypeError(f"Expected swanlab.Callback, got {type(cb).__name__}")

            if cb.name in self._callbacks:
                console.warning(f"Callback '{cb.name}' is already registered and will be overwritten.")

            self._callbacks[cb.name] = cb

    def remove_callback(self, name: str) -> None:
        """根据名称移除指定的回调函数"""
        if name in self._callbacks:
            del self._callbacks[name]

    @property
    def registered_callbacks(self) -> List[Callback]:
        """返回当前所有已注册的回调列表"""
        return list(self._callbacks.values())

    # =========================================================================
    # 动态事件分发引擎 (Dynamic Event Dispatcher)
    # 拦截所有对未定义方法的调用，自动转发给注册的 callbacks
    # =========================================================================
    def __getattr__(self, name: str):
        # 1. 安全检查：只有基类 Callback 中真实定义的方法才允许被动态调用
        if hasattr(Callback, name) and callable(getattr(Callback, name)):
            # 2. 动态生成带容错机制的代理函数
            def dispatcher(*args, **kwargs):
                for cb in self._callbacks.values():
                    # 从具体的回调实例中获取对应名称的方法并执行
                    with safe.block(message=f"Error executing '{name}' in callback '{cb.name}'"):
                        getattr(cb, name)(*args, **kwargs)

            # 3. 性能优化：将生成好的代理函数缓存到当前实例上
            # 这样第二次调用 callbacker.on_log 时，就变成了普通的 O(1) 属性访问，性能极高！
            setattr(self, name, dispatcher)
            return dispatcher

        # 如果调用的既不是已定义方法，也不是基类中的回调方法，则抛出标准异常
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")


# 全局 CallBacker，主要用处是存储全局注册的回调函数
global_callbacker = _CallbackManager()

# =========================================================================
# IDE 类型欺骗 (Type Hinting Magic)
# 这一段代码在真实运行时根本不会被执行，它完全是写给 PyCharm/VSCode 看的
# =========================================================================
if TYPE_CHECKING:
    # 让 IDE 认为全局的 callbacker 同时拥有 CallbackManager 和 Callback 的所有方法
    class CallbackManager(_CallbackManager, Callback):
        @property
        def name(self) -> str:
            return "FakeCallback"

else:
    CallbackManager = _CallbackManager


def create_callback_manager(callbacks: Optional[CallbacksType]) -> CallbackManager:
    """
    创建一个新的回调管理器，继承全局回调的同时，支持注入局部回调（当前run有效的回调）
    优先级是局部回调大于全局回调
    """
    cm = CallbackManager()
    cm.merge_callbacks(global_callbacker.registered_callbacks)
    cm.merge_callbacks(callbacks)
    return cm


__all__ = ["CallbackManager", "global_callbacker", "create_callback_manager"]
