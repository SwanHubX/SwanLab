#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/12 14:20
@File: operator.py
@IDE: pycharm
@Description:
    回调函数操作员，批量处理回调函数的调用
"""
from typing import List, Union, Callable
from .callback import SwanLabRunCallback
from ..settings import SwanDataSettings


def fmt_return(func):
    """
    对返回值进行格式化，如果返回值是一个列表，那么就返回最后一个不是None的值，否则直接返回
    如果所有的值都是None，那么返回None
    """

    def wrapper(*args, **kwargs):
        rets = func(*args, **kwargs)
        if isinstance(rets, list):
            for ret in rets:
                if ret is not None:
                    return ret
            return None
        return rets

    return wrapper


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
        self.callbacks: List[SwanLabRunCallback] = callbacks if callbacks is not None else []

    def inject(self, settings: SwanDataSettings):
        self.settings = settings
        return [callback.inject(settings) for callback in self.callbacks]

    def add_callback(self, callback: SwanLabRunCallback):
        self.callbacks.append(callback)

    @fmt_return
    def before_init_experiment(self, *args, **kwargs):
        return [callback.before_init_experiment(*args, **kwargs) for callback in self.callbacks]

    @fmt_return
    def before_init_project(self, *args, **kwargs):
        return [callback.before_init_project(*args, **kwargs) for callback in self.callbacks]

    @fmt_return
    def on_log(self, *args, **kwargs):
        return [callback.on_log(*args, **kwargs) for callback in self.callbacks]

    @fmt_return
    def on_train_begin(self):
        return [callback.on_train_begin() for callback in self.callbacks]

    @fmt_return
    def on_train_end(self, *args, **kwargs):
        return [callback.on_train_end(*args, **kwargs) for callback in self.callbacks]

    @fmt_return
    def on_metric_create(self, *args, **kwargs):
        return [callback.on_metric_create(*args, **kwargs) for callback in self.callbacks]

    @fmt_return
    def on_column_create(self, *args, **kwargs):
        return [callback.on_column_create(*args, **kwargs) for callback in self.callbacks]
