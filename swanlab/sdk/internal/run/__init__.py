"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:29
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理
2. 运行、实验上下文管理
3. 触发回调
"""

from functools import cached_property, singledispatchmethod
from typing import Any, Dict, Optional, Union, get_args

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType

from .callbackers import CloudCallback, LocalCallback, OfflineCallback
from .data.transforms import Text

__all__ = ["SwanLabRun", "CloudCallback", "LocalCallback", "OfflineCallback"]


class SwanLabRun:
    def __init__(self, ctx: "RunContext"):
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
        # 2. 展开字典，例如 {"a": {"b": {"c": 1}}} -> {"a/b/c": 1}
        data = flatten_dict(data)
        # 3. 依次处理
        for key, value in data.items():
            self._dispatch_log(value=value, key=key, step=step)

    @singledispatchmethod
    def _dispatch_log(self, value: Any, key: str, step: Optional[int] = None):
        """默认的回退处理逻辑
        此部分用于处理基础标量数据类型
        """
        # self.log_scalar(key=key, step=step, value=value)
        pass

    @_dispatch_log.register
    def _(self, value: Text, key: str, step: Optional[int]):
        """Log一个Text"""
        return self.log_text(key=key, step=step, data=value)

    def log_text(self, key: str, data: Union[str, Text], caption: Optional[str] = None, step: Optional[int] = None):
        """
        Log a text record to the current run.
        :param key: The key of the record.
        :param data: The text content or a Text object.
        :param caption: An optional caption for the text.
        :param step: The global step at which the data was logged. Can be None if not explicitly tracked.
        """
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


def flatten_dict(
    d: Dict[str, Any], parent_key: str = "", parent_dict: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    展开字典，例如 {"a": {"b": {"c": 1}}} -> {"a/b/c": 1}
    如果出现重复键名，根据顺序，顺序靠后的键名会覆盖靠前的键名
    :param d: 待展开的字典
    :param parent_key: 父级键名，用于构建新键名
    :param parent_dict: 父级字典，用于存储展开后的结果
    :return: 展开后的字典
    """
    # 顶层调用时初始化字典（避免可变默认参数陷阱）
    if parent_dict is None:
        parent_dict = {}

    for k, v in d.items():
        new_key = f"{parent_key}/{k}" if parent_key else k

        if isinstance(v, dict):
            # 递归调用，将同一个 parent_dict 引用传递下去
            flatten_dict(v, new_key, parent_dict)
        else:
            # 检查冲突并警告
            if new_key in parent_dict:
                console.warning(f"Duplicate key found: '{new_key}'. The latter value will overwrite the former one.")
            # 直接在共享的字典上赋值
            parent_dict[new_key] = v

    return parent_dict
