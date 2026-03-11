"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 17:12
@description: SwanLab 运行时回调函数
"""

from typing import TYPE_CHECKING, Dict, Iterable, List

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.utils.callbacker import SwanLabCallback


class _CallbackManager:
    """
    SwanLab 回调函数的管理器。
    维护当前上下文中的全局回调状态，支持批量合并与指定移除。
    """

    def __init__(self):
        # 使用字典保证顺序并天然去重
        self._callbacks: Dict[str, SwanLabCallback] = {}

    def merge_callbacks(self, callbacks: Iterable[SwanLabCallback]) -> None:
        """批量合并回调函数到当前管理器中"""
        if not callbacks:
            return

        for cb in callbacks:
            if not isinstance(cb, SwanLabCallback):
                raise TypeError(f"Expected SwanLabCallback, got {type(cb).__name__}")

            if cb.name in self._callbacks:
                console.warning(f"Callback '{cb.name}' is already registered and will be overwritten.")

            self._callbacks[cb.name] = cb

    def remove_callback(self, name: str) -> None:
        """根据名称移除指定的回调函数"""
        if name in self._callbacks:
            del self._callbacks[name]

    @property
    def registered_callbacks(self) -> List[SwanLabCallback]:
        """返回当前所有已注册的回调列表"""
        return list(self._callbacks.values())

    # =========================================================================
    # 动态事件分发引擎 (Dynamic Event Dispatcher)
    # 拦截所有对未定义方法的调用，自动转发给注册的 callbacks
    # =========================================================================
    def __getattr__(self, name: str):
        # 1. 安全检查：只有基类 SwanLabCallback 中真实定义的方法才允许被动态调用
        if hasattr(SwanLabCallback, name) and callable(getattr(SwanLabCallback, name)):
            # 2. 动态生成带容错机制的代理函数
            def dispatcher(*args, **kwargs):
                for cb in self._callbacks.values():
                    try:
                        # 从具体的回调实例中获取对应名称的方法并执行
                        getattr(cb, name)(*args, **kwargs)
                    except Exception as e:
                        console.error(f"Error executing '{name}' in callback '{cb.name}': {e}")

            # 3. 性能优化：将生成好的代理函数缓存到当前实例上
            # 这样第二次调用 callbacker.on_log 时，就变成了普通的 O(1) 属性访问，性能极高！
            setattr(self, name, dispatcher)
            return dispatcher

        # 如果调用的既不是已定义方法，也不是基类中的回调方法，则抛出标准异常
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")


# =========================================================================
# IDE 类型欺骗 (Type Hinting Magic)
# 这一段代码在真实运行时根本不会被执行，它完全是写给 PyCharm/VSCode 看的
# =========================================================================
if TYPE_CHECKING:
    # 让 IDE 认为全局的 callbacker 同时拥有 CallbackManager 和 SwanLabCallback 的所有方法
    class CallbackManager(_CallbackManager, SwanLabCallback):
        @property
        def name(self) -> str:
            return "FakeCallback"

    callbacker: CallbackManager
else:
    CallbackManager = _CallbackManager
    # 全局 CallBacker
    callbacker = CallbackManager()

__all__ = ["callbacker", "CallbackManager"]
