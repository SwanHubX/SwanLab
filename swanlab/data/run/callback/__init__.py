#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/5 16:27
@File: __init__.py
@IDE: pycharm
@Description:
    回调函数注册抽象模块
"""
from typing import Callable
from abc import ABC, abstractmethod
from .utils import U
from .models import *
import atexit
import sys


class SwanLabRunCallback(ABC, U):
    """
    SwanLabRunCallback，回调函数注册类，所有以`on_`和`before_`开头的函数都会在对应的时机被调用
    为了方便管理：
    1. `_`开头的函数为内部函数，不会被调用，且写在最开头
    2. 所有回调按照逻辑上的触发顺序排列
    3. 带有from_*后缀的回调函数代表调用者来自其他地方，比如config、operator等，这将通过settings对象传递
    4. 所有回调不要求全部实现，只需实现需要的回调即可
    """

    def _register_sys_callback(self):
        """
        注册系统回调，内部使用
        """
        sys.excepthook = self._except_handler
        atexit.register(self._clean_handler)

    def _unregister_sys_callback(self):
        """
        注销系统回调，内部使用
        """
        sys.excepthook = sys.__excepthook__
        atexit.unregister(self._clean_handler)

    def _clean_handler(self):
        """
        正常退出清理函数，此函数调用`run.finish`
        """
        pass

    def _except_handler(self, tp, val, tb):
        """
        异常退出清理函数
        """
        pass

    def on_init(self, proj_name: str, workspace: str, logdir: str = None):
        """
        执行`swanlab.init`时调用
        此时运行时环境变量没有被设置，此时修改环境变量还是有效的
        :param logdir: str, 用户设置的日志目录
        :param proj_name: str, 项目名称
        :param workspace: str, 工作空间
        """
        pass

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        suffix: str,
        setter: Callable[[str, str, str, str], None],
    ):
        """
        在初始化实验之前调用，此时SwanLabRun已经初始化完毕
        :param run_id: str, SwanLabRun的运行id
        :param exp_name: str, 实验名称
        :param description: str, 实验描述
        :param num: int, 历史实验数量
        :param suffix: str, 实验后缀
        :param setter: Callable[[str, str, str, str], None], 设置实验信息的函数，在这里设置实验信息
        """
        pass

    def on_run(self):
        """
        SwanLabRun初始化完毕时调用
        """
        pass

    def on_run_error_from_operator(self, e: OperateErrorInfo):
        """
        SwanLabRun初始化错误时被操作员调用
        """
        pass

    def on_runtime_info_update(self, r: RuntimeInfo):
        """
        运行时信息更新时调用
        :param r: RuntimeInfo, 运行时信息
        """
        pass

    def on_log(self):
        """
        每次执行swanlab.log时调用
        """
        pass

    def on_column_create(self, column_info: ColumnInfo):
        """
        列创建回调函数,新增列信息时调用
        """
        pass

    def on_metric_create(self, metric_info: MetricInfo):
        """
        指标创建回调函数,新增指标信息时调用
        """
        pass

    def on_stop(self, error: str = None):
        """
        训练结束时的回调函数
        """
        pass

    @abstractmethod
    def __str__(self) -> str:
        """
        返回当前回调函数的名称，这条应该是一个全局唯一的标识
        在operator中会用到这个名称，必须唯一
        """
        pass
