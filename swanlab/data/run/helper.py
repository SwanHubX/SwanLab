#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/12 14:20
@File: helper.py
@IDE: pycharm
@Description:
    回调函数操作员，批量处理回调函数的调用
"""
from enum import Enum
from typing import List, Union, Dict, Any, Tuple

from swanlab.data.run.webhook import try_send_webhook
from swanlab.toolkit import SwanKitCallback, MetricInfo, ColumnInfo, RuntimeInfo

OperatorReturnType = Dict[str, Any]


class SwanLabRunOperator(SwanKitCallback):
    """
    The SwanLabRunOperator is used to batch process callback instances, triggering them respectively in the order of
    injection when each occasion occurs. In design, we aim for isolation between each instance. However, if there are
    cases of shared scope, parameters, or variables between instances, developers should handle them appropriately,
    as SwanLabRunOperator itself will not address these situations for now.
    """

    def __init__(self, callbacks: Union[SwanKitCallback, List[SwanKitCallback]] = None):
        super(SwanLabRunOperator, self).__init__()
        callbacks = [callbacks] if isinstance(callbacks, SwanKitCallback) else callbacks
        self.callbacks = {}
        if callbacks is not None:
            for callback in callbacks:
                self.add_callback(callback)

    @property
    def disabled(self):
        """
        判断是否所有回调函数都已经被禁用
        FIXME 实际上这里与settings的一些属性有关，目前没问题，但是未来可能要改
        """
        return len(self.callbacks) == 0

    def add_callback(self, callback: SwanKitCallback):
        if not isinstance(callback, SwanKitCallback):
            raise TypeError(f"Unsupported callback type: {type(callback)}")
        if str(callback) == str(self) or callback in self.callbacks:
            raise ValueError(f"Cannot add the same callback instance: {callback}")
        self.callbacks[str(callback)] = callback

    def __str__(self):
        return "SwanLabRunOperator"

    def __run_all(self, method: str, *args, **kwargs):
        return {name: getattr(callback, method)(*args, **kwargs) for name, callback in self.callbacks.items()}

    @classmethod
    def parse_return(cls, ret: OperatorReturnType, key: str = None):
        """
        解析返回值，选择不为None的返回值，如果都为None，则返回None
        如果key不为None，则返回对应key的返回值
        :param ret: 返回值
        :param key: 返回值的key
        :return:
        """
        if key is not None:
            return ret[key]
        return next((v for v in ret.values() if v is not None), None)

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs) -> OperatorReturnType:
        return self.__run_all(
            "on_init",
            proj_name,
            workspace,
            logdir=logdir,
            *args,
            **kwargs,
        )

    def before_run(self, *args, **kwargs):
        return self.__run_all("before_run", *args, **kwargs)

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        return self.__run_all(
            "before_init_experiment",
            run_id,
            exp_name,
            description,
            colors,
            *args,
            **kwargs,
        )

    def on_run(self, *args, **kwargs):
        self.__run_all("on_run", *args, **kwargs)
        try_send_webhook()

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        return self.__run_all("on_runtime_info_update", r, *args, **kwargs)

    def on_log(self, *args, **kwargs):
        return self.__run_all("on_log", *args, **kwargs)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        return self.__run_all("on_metric_create", metric_info, *args, **kwargs)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        return self.__run_all("on_column_create", column_info, *args, **kwargs)

    def on_stop(self, error: str = None, epoch: int = None, *args, **kwargs):
        r = self.__run_all("on_stop", error=error, epoch=epoch, *args, **kwargs)
        # 清空所有注册的回调函数
        self.callbacks.clear()
        return r


class SwanLabRunState(Enum):
    """SwanLabRunState is an enumeration class that represents the state of the experiment.
    We Recommend that you use this enumeration class to represent the state of the experiment.
    """

    NOT_STARTED = -2
    SUCCESS = 1
    CRASHED = -1
    RUNNING = 0


def monitor_interval(count: int) -> float:
    """
    根据传入的参数（本质是调用次数）选择间隔
    :param count: (调用)次数
    """
    if count < 10:
        return 10
    elif count < 50:
        return 30
    else:
        return 60
