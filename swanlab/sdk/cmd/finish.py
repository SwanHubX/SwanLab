"""
@author: cunyue
@file: finish.py
@time: 2026/3/6 21:49
@description: SwanLab SDK 结束当前运行
"""

import atexit
import sys
import traceback
from types import TracebackType
from typing import Optional, Type

from swanlab.sdk.cmd.helper import with_cmd_lock
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import get_run, has_run
from swanlab.sdk.typings.run import FinishType


@with_cmd_lock
def finish(state: FinishType = "success", error: Optional[str] = None):
    """
    Finish the current run and close the current experiment
    Normally, swanlab will run this function automatically,
    but you can also execute it manually and mark the experiment as 'completed'.
    Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
    If you mark the experiment as 'CRASHED' manually, `error` must be provided.
    """
    if not has_run():
        console.error("No active SwanLabRun. Call swanlab.init() first.")
        return
    run = get_run()
    run.finish(state, error)


def atexit_finish():
    """
    全局退出时自动结束当前运行，此时代码正常执行完毕
    """
    if not has_run():
        return
    console.debug("SwanLab Run is finishing at exit...")
    run = get_run()
    run.finish()


atexit.register(atexit_finish)


def swanlab_excepthook(tp: Type[BaseException], val: BaseException, tb: Optional[TracebackType]):
    """全局异常捕获，用于将实验标记为 crashed"""
    try:
        if not has_run():
            return
        state: FinishType = "crashed"
        if tp is KeyboardInterrupt:
            console.info("KeyboardInterrupt by user")
            state = "aborted"
        else:
            console.info("Error happened while training")
        # 生成错误堆栈
        full_error_msg = "".join(traceback.format_exception(tp, val, tb))

        # 打印错误堆栈
        run = get_run()
        run.finish(state=state, error=full_error_msg)

    except Exception as e:
        console.error(f"SwanLab  failed to handle excepthook: {e}", file=sys.stderr)
    finally:
        # _original_excepthook = sys.excepthook
        # 不要用动态保存的 _original_excepthook，直接调用 Python 底层 C 实现的 sys.__excepthook__
        # 确保异常信息被正确打印
        sys.__excepthook__(tp, val, tb)


sys.excepthook = swanlab_excepthook
