#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/12 14:20
@File: operator.py
@IDE: pycharm
@Description:
    回调函数操作员，批量处理回调函数的调用
"""
from typing import List, Union
from .callback import SwanLabRunCallback
from ..settings import SwanDataSettings


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
        [callback.inject(settings) for callback in self.callbacks]

    def add_callback(self, callback: SwanLabRunCallback):
        self.callbacks.append(callback)

    def before_init_project(self, *args, **kwargs):
        [callback.before_init_project(*args, **kwargs) for callback in self.callbacks]

    def on_train_begin(self, *args, **kwargs):
        [callback.on_train_begin(*args, **kwargs) for callback in self.callbacks]

    def on_train_end(self, *args, **kwargs):
        [callback.on_train_end(*args, **kwargs) for callback in self.callbacks]

    def on_metric_create(self, *args, **kwargs):
        [callback.on_metric_create(*args, **kwargs) for callback in self.callbacks]

    def on_column_create(self, *args, **kwargs):
        [callback.on_column_create(*args, **kwargs) for callback in self.callbacks]
