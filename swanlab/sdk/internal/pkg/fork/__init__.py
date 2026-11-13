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
- 生命周期安全：bound method 回调以弱引用存储，注册不会延长实例生命周期；
  回调内禁止获取实例级锁（fork 时该锁可能由不会存在于子进程的线程持有）
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

import inspect
import os
import threading
import weakref
from typing import Callable, List, Optional

# 回调注册表：fork 后在子进程中执行的清理函数。
# bound method 以 WeakMethod 存储：注册不得延长实例生命周期，
# 否则实例上的 weakref.finalize 等最后保险永远不会触发。
# 模块级函数与 lambda 保持强引用（其生命周期本就不随实例结束）。
# 使用锁保护，因为 register() 可能在多线程环境中被调用
_callbacks: List[object] = []
_lock = threading.Lock()


def _make_entry(callback: Callable[[], None]) -> object:
    if inspect.ismethod(callback):
        return weakref.WeakMethod(callback)
    return callback


def _resolve_entry(entry: object) -> Optional[Callable[[], None]]:
    if isinstance(entry, weakref.WeakMethod):
        return entry()
    return entry  # type: ignore[return-value]


def _prune_dead() -> None:
    """清除已死亡的弱引用条目。需持锁调用。"""
    _callbacks[:] = [e for e in _callbacks if _resolve_entry(e) is not None]


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
    特别注意：回调内不得获取可能被其他线程持有、而该线程不会存在于
    子进程中的锁，否则子进程会在 fork 返回前永久阻塞。

    bound method 以弱引用存储：实例被回收后回调自动失效；注册方仍应
    在生命周期结束时调用 ``unregister`` 即时移除条目。

    :param callback: 无参回调函数，在 fork 后的子进程中执行
    """
    with _lock:
        _prune_dead()
        _callbacks.append(_make_entry(callback))


def unregister(callback: Callable[[], None]) -> None:
    """注销已注册的 fork 回调。

    用于 Run 等对象在 finish 时清理自己注册的回调，避免回调列表持续增长。
    如果 callback 未注册过则静默忽略。

    :param callback: 之前通过 ``register()`` 注册的回调函数
    """
    with _lock:
        _prune_dead()
        for index, entry in enumerate(_callbacks):
            if _resolve_entry(entry) == callback:
                del _callbacks[index]
                return


def _before_fork() -> None:
    """``os.register_at_fork(before=...)`` 的处理器。

    在 fork 之前获取锁，确保没有其他线程正在修改 _callbacks。
    """
    _lock.acquire()


def _after_in_parent() -> None:
    """``os.register_at_fork(after_in_parent=...)`` 的处理器。

    fork 完成后在父进程中释放锁。
    """
    _lock.release()


def _after_in_child() -> None:
    """``os.register_at_fork(after_in_child=...)`` 的处理器。

    在 fork 后的子进程中执行所有已注册的回调，然后重置锁。
    回调执行顺序与注册顺序一致（FIFO）；弱引用已死亡的条目被跳过
    （父进程侧的注册表不因 fork 而改变）。

    因为 _before_fork 已获取锁，此处锁处于持有状态。
    直接重置为新锁，避免继承父进程的锁状态。
    """
    # _before_fork 已获取锁，此时 _callbacks 一定处于一致状态
    try:
        for entry in list(_callbacks):
            cb = _resolve_entry(entry)
            if cb is not None:
                cb()
    finally:
        global _lock
        _lock = threading.Lock()


# 注册 os.register_at_fork 回调（仅 POSIX 平台可用）
if hasattr(os, "register_at_fork"):
    os.register_at_fork(before=_before_fork, after_in_parent=_after_in_parent, after_in_child=_after_in_child)
