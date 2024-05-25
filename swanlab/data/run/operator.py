#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/12 14:20
@File: operator.py
@IDE: pycharm
@Description:
    回调函数操作员，批量处理回调函数的调用
"""
from typing import List, Union, Dict, Any, Callable
from .callback import SwanLabRunCallback, MetricInfo, ColumnInfo
from ..settings import SwanDataSettings

OperatorReturnType = Dict[str, Any]


class SwanLabRunOperator(SwanLabRunCallback):
    """
    The SwanLabRunOperator is used to batch process callback instances, triggering them respectively in the order of
    injection when each occasion occurs. In design, we aim for isolation between each instance. However, if there are
    cases of shared scope, parameters, or variables between instances, developers should handle them appropriately,
    as SwanLabRunOperator itself will not address these situations for now.
    """

    def __init__(self, callbacks: Union[SwanLabRunCallback, List[SwanLabRunCallback]] = None):
        super(SwanLabRunOperator, self).__init__()
        callbacks = [callbacks] if isinstance(callbacks, SwanLabRunCallback) else callbacks
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

    def add_callback(self, callback: SwanLabRunCallback):
        if not isinstance(callback, SwanLabRunCallback):
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
        解析返回值，选择不为None的返回值
        如果key不为None，则返回对应key的返回值
        :param ret: 返回值
        :param key: 返回值的key
        :return:
        """
        if key is not None:
            return ret[key]
        return next((v for v in ret.values() if v is not None), None)

    def on_init(self, proj_name: str, workspace: str, logdir: str = None) -> OperatorReturnType:
        return self.__run_all("on_init", proj_name, workspace, logdir=logdir)

    def inject(self, settings: SwanDataSettings) -> OperatorReturnType:
        self.settings = settings
        return {name: callback.inject(settings) for name, callback in self.callbacks.items()}

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        suffix: str,
        setter: Callable[[str, str, str, str], None]
    ):
        return self.__run_all(
            "before_init_experiment",
            run_id,
            exp_name,
            description,
            num,
            suffix,
            setter
        )

    def on_run(self):
        return self.__run_all("on_run")

    def on_log(self):
        return self.__run_all("on_log")

    def on_metric_create(self, metric_info: MetricInfo):
        return self.__run_all("on_metric_create", metric_info)

    def on_column_create(self, column_info: ColumnInfo):
        return self.__run_all("on_column_create", column_info)

    def on_stop(self, error: str = None):
        return self.__run_all("on_stop", error)
