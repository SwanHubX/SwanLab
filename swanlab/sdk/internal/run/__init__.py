"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:29
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理
2. 运行、实验上下文管理
3. 触发回调
"""

from functools import cached_property
from typing import Any, Dict, Optional, get_args

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType


class SwanLabRun:
    def __init__(self, ctx: RunContext):
        self._ctx = ctx

    @cached_property
    def id(self) -> str:
        """
        This property returns the unique identifier of the current SwanLab run.
        """
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    def log(self, data: Dict[str, Any], step: Optional[int] = None):
        """
        Log a row of data to the current run. Unlike `swanlab.log`, this api will be called directly based on the
        SwanRun instance, removing the initialization process. Of course, after you call the finish method,
        this method will be banned from calling.

        :param data: The data to log.
        :param step: The global step at which the data was logged. Can be None if not explicitly tracked.
        """
        # 1. 类型检查
        if not isinstance(data, dict):
            console.error("Log data must be a dict, but got {}. SwanLab will ignore records it.".format(type(data)))
            return
        ...

    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """
        Finish the current run and close the current experiment
        Normally, swanlab will run this function automatically,
        but you can also execute it manually and mark the experiment as 'completed'.
        Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
        After calling this function, you can re-run `swanlab.init` to start a new experiment.

        :param state: The state of the experiment, it can be 'success', 'crashed' or 'aborted'.
        :param error: The error message when the experiment is marked as 'crashed'. If not 'crashed', it should be None.
        """
        # 类型检查
        state = state.lower()  # type: ignore
        if state not in get_args(FinishType):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return
        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"
        ...
