"""
@author: cunyue
@file: __init__.py
@time: 2026/4/19
@description: SwanLab SDK fork 感知模块，提供：

1. ``current_pid()`` — 获取当前进程 PID
2. ``is_forked(pre_pid)`` — 判断当前进程是否为 fork 子进程
3. ``register(callback)`` — 注册 fork 后子进程的清理回调
4. ``unregister(callback)`` — 注销已注册的回调
5. ``_after_in_child()`` — ``os.register_at_fork(after_in_child=...)`` 的处理器

设计原则：

- 回调注册表：各组件通过 ``register()`` 注册自己的 fork 清理逻辑，而非内联 ``register_at_fork`` 调用
- 显式优于隐式：fork 后的清理通过回调显式触发，而非依赖属性计算的副作用
- 面向未来：当 ``swanlab-core`` 上线后，可注册重连回调替代简单的 ``clear_run()``

使用示例::

    from swanlab.sdk.internal.pkg import fork

    # 获取当前 PID（用于记录，后续传给 is_forked）
    pid = fork.current_pid()

    # 判断是否 fork
    if fork.is_forked(pid):
        ...

    # 注册 fork 回调
    fork.register(lambda: some_cleanup())

    # 注销回调
    fork.unregister(some_cleanup)
"""

import os
import threading
from typing import Callable, List

# 回调列表：fork 后在子进程中执行的清理函数
# 使用锁保护，因为 register() 可能在多线程环境中被调用
_callbacks: List[Callable[[], None]] = []
_lock = threading.Lock()


def current_pid() -> int:
    """获取当前进程 PID。

    所有需要记录 PID 的地方应通过此函数获取，而非直接调用 ``os.getpid()``，
    确保 fork 模块为 PID 操作的唯一入口。

    :return: 当前进程 PID
    """
    return os.getpid()


def is_forked(pre_pid: int) -> bool:
    """判断当前进程是否为 fork 子进程。

    通过比较当前 PID 与传入的 pre_pid 来判断。
    调用方需在对象创建时通过 ``current_pid()`` 记录 PID，后续传给此函数检测。

    :param pre_pid: 需要比较的 PID，通常为对象创建时的 ``current_pid()`` 返回值
    :return: 如果当前进程是 fork 子进程则返回 True
    """
    return os.getpid() != pre_pid


def register(callback: Callable[[], None]) -> None:
    """注册 fork 后子进程的清理回调。

    回调在 ``os.register_at_fork(after_in_child=...)`` 触发时执行，
    即 fork 后子进程的第一时间。回调应尽量轻量，避免阻塞或抛出异常。

    :param callback: 无参回调函数，在 fork 后的子进程中执行
    """
    with _lock:
        _callbacks.append(callback)


def unregister(callback: Callable[[], None]) -> None:
    """注销已注册的 fork 回调。

    用于 Run 等对象在 finish 时清理自己注册的回调，避免回调列表持续增长。
    如果 callback 未注册过则静默忽略。

    :param callback: 之前通过 ``register()`` 注册的回调函数
    """
    with _lock:
        try:
            _callbacks.remove(callback)
        except ValueError:
            pass


def _after_in_child() -> None:
    """``os.register_at_fork(after_in_child=...)`` 的处理器。

    在 fork 后的子进程中执行所有已注册的回调，然后替换锁。
    回调执行顺序与注册顺序一致（FIFO）。

    替换锁的原因：fork 时父进程的 _lock 可能被其他线程持有，
    子进程继承的旧锁永远无法释放，后续 register/unregister 会死锁。
    """
    # 子进程中只有调用 fork 的线程存在，此时读 _callbacks 是安全的
    for cb in _callbacks:
        cb()
    global _lock
    _lock = threading.Lock()


# 注册 os.register_at_fork 回调（仅 POSIX 平台可用）
if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_in_child)
