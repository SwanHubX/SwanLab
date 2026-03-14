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

from swanlab.sdk.cmd.helper import with_cmd_lock, with_run
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import get_run, has_run
from swanlab.sdk.typings.run import FinishType


@with_cmd_lock
@with_run("finish")
def finish(state: FinishType = "success", error: Optional[str] = None):
    """Finish the current run and close the experiment.

    This function safely closes the current run and waits for all logs to be flushed.
    SwanLab automatically calls this function at program exit, but you can call it
    manually to mark the experiment as completed with a specific state.

    :param state: Final state of the run. Must be one of: "success", "crashed", "aborted".
        Defaults to "success".

    :param error: Error message if state is "crashed". Required when state="crashed".

    :raises RuntimeError: If called without an active run.

    Examples:

        Finish a successful run:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Mark run as crashed with error message:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> try:
        ...     # training code
        ...     raise ValueError("Training failed")
        ... except Exception as e:
        ...     swanlab.finish(state="crashed", error=str(e))

        Mark run as aborted:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish(state="aborted")
    """
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
        console.error(f"SwanLab failed to handle excepthook: {e}", file=sys.stderr)
    finally:
        # _original_excepthook = sys.excepthook
        # 不要用动态保存的 _original_excepthook，直接调用 Python 底层 C 实现的 sys.__excepthook__
        # 确保异常信息被正确打印
        sys.__excepthook__(tp, val, tb)


sys.excepthook = swanlab_excepthook
