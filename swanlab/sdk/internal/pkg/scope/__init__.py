"""
@author: cunyue
@file: __init__.py
@time: 2026/3/8 13:55
@description: 有些时候，当我们需要获取到一个深层函数的返回值或者临时信息
为了从底层函数捞出一个临时信息，而被迫修改中间所有链路的返回值或参数签名（通常被称为“透传污染”），绝对是一种极具破坏性的反模式
在 Python 中，解决这个问题的现代且线程/协程安全的方案，是结合 contextvars（上下文变量）和魔法方法（__enter__, __exit__）来自己实现一个上下文收集器对象。
这个对象被称为 Scope，它允许我们在不修改函数签名的情况下，从深层函数中获取返回值或临时信息。
于此同时，Scope支持嵌套，允许我们从深层函数中获取多个临时信息，或者在嵌套调用中获取不同层次的临时信息。
"""

import contextvars
from typing import Any, Dict, Optional

# 核心：定义全局的 ContextVar，存储当前处于激活状态的 Scope 实例
# 它在多线程和 asyncio 协程下是绝对隔离且安全的
_scope_ctx: contextvars.ContextVar[Optional["Scope"]] = contextvars.ContextVar("swanlab_scope", default=None)


class Scope:
    """
    上下文作用域管理器。
    """

    def __init__(self, bubble_up: bool = False):
        """
        初始化一个独立的作用域。

        :param bubble_up: 是否在退出当前 Scope 时，将收集到的数据合并（冒泡）到父级 Scope 中。
                          如果你希望顶层 Scope 能在最后拿到所有子 Scope 收集的汇总数据，请设为 True。
        """
        self.data: Dict[str, Any] = {}
        self.bubble_up = bubble_up
        self._token: Optional[contextvars.Token] = None
        self._parent: Optional["Scope"] = None

    def __enter__(self) -> "Scope":
        # 1. 记录可能存在的父级 Scope（用于嵌套链式查找或冒泡合并）
        self._parent = _scope_ctx.get()

        # 2. 夺取控制权，将当前实例设置为全局激活的上下文
        self._token = _scope_ctx.set(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # 1. 退出时，安全地归还上下文控制权给父级
        if self._token is not None:
            _scope_ctx.reset(self._token)
            self._token = None

        # 2. 退出清理：如果开启了冒泡，并且存在父级 Scope，则把自己的数据上报合并给父级
        if self.bubble_up and self._parent is not None:
            self._parent.data.update(self.data)

    def set(self, key: str, value: Any) -> None:
        """向当前 Scope 中写入数据"""
        self.data[key] = value

    def get(self, key: str, default: Any = None) -> Any:
        """
        从当前 Scope 中读取数据。
        【嵌套特性】：如果当前层级没有这个数据，会自动向上级 Scope 追溯寻找。
        """
        if key in self.data:
            return self.data[key]
        if self._parent is not None:
            return self._parent.get(key, default)
        return default


# ==============================================================================
# 暴露给业务代码的快捷函数，避免业务代码直接操作全局变量 _scope_ctx
# ==============================================================================


def set_context(key: str, value: Any) -> None:
    """
    【隔空写入】跨层级写入上下文数据。
    如果当前不在任何 Scope 的 with 块中，则静默忽略（非常安全）。
    """
    current_scope = _scope_ctx.get()
    if current_scope is not None:
        current_scope.set(key, value)


def get_context(key: str, default: Any = None) -> Any:
    """
    【隔空读取】跨层级读取上下文数据。
    如果当前不在任何 Scope 的 with 块中，则返回 default。
    """
    current_scope = _scope_ctx.get()
    if current_scope is not None:
        return current_scope.get(key, default)
    return default
