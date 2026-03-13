"""
@author: cunyue
@file: finish.py
@time: 2026/3/6 21:49
@description: SwanLab SDK 结束当前运行
"""

from typing import Optional

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
