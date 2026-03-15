"""
@author: cunyue
@file: run.py
@time: 2026/3/14
@description: SwanLab SDK run methods

_make_run_cmd 将 SwanLabRun 上的方法包装为顶层 cmd 函数：
  - __wrapped__ 指向 SwanLabRun 原方法，help() / inspect 可追溯
  - __doc__ 复用原方法文档，无需重复维护

新增公开方法时，在 SwanLabRun 上实现后，在此文件末尾追加一行：
    xxx = _make_run_cmd("xxx")
并同步更新 swanlab/__init__.pyi 中的函数声明。
"""

from typing import Any, Callable

from swanlab.sdk.cmd.helper import with_cmd_lock, with_run
from swanlab.sdk.internal.run import SwanLabRun, get_run


def _make_run_cmd(method_name: str) -> Callable:
    run_method = getattr(SwanLabRun, method_name)

    @with_cmd_lock
    @with_run(method_name)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        return getattr(get_run(), method_name)(*args, **kwargs)

    wrapper.__name__ = method_name
    wrapper.__wrapped__ = run_method  # type: ignore[attr-defined]
    wrapper.__doc__ = run_method.__doc__
    return wrapper


# ── 每新增一个公开方法，在此追加一行 ──────────────────────────────────────────
log = _make_run_cmd("log")
log_text = _make_run_cmd("log_text")
log_image = _make_run_cmd("log_image")
log_audio = _make_run_cmd("log_audio")
log_video = _make_run_cmd("log_video")
define_scalar = _make_run_cmd("define_scalar")
finish = _make_run_cmd("finish")
