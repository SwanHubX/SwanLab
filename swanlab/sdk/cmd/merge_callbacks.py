"""
@author: cunyue
@file: merge_callbacks.py
@time: 2026/4/30
@description: swanlab.merge_callbacks 方法，合并自定义回调函数
"""

from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.context.components.callbacker import global_callbacker
from swanlab.sdk.typings.context import CallbacksType

__all__ = ["merge_callbacks"]


@with_cmd_lock
@without_run("merge_callbacks")
def merge_callbacks(callbacks: CallbacksType) -> None:
    """Merge custom callbacks into the global SwanLab callback registry.

    This function allows you to register callbacks before initializing a run.
    It must be called before `swanlab.init()`.

    :param callbacks: Custom callbacks to merge. Can be a single Callback object or an iterable of Callback objects.

    Examples:

        Merge a single callback:

        >>> import swanlab
        >>> from swanlab import Callback
        >>> class MyCallback(Callback):
        ...     @property
        ...     def name(self) -> str:
        ...         return "my_callback"
        ...     def on_run_initialized(self, run_dir, path):
        ...         print("Run initialized!")
        >>> swanlab.merge_callbacks(MyCallback())
        >>> swanlab.init()

        Merge multiple callbacks:

        >>> class AnotherCallback(Callback):
        ...     @property
        ...     def name(self) -> str:
        ...         return "another_callback"
        ...     def on_run_initialized(self, run_dir, path):
        ...         print("Run initialized by another callback!")
        >>> swanlab.merge_callbacks([MyCallback(), AnotherCallback()])
        >>> swanlab.init()
    """
    global_callbacker.merge_callbacks(callbacks)
