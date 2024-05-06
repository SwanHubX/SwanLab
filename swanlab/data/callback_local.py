#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback_local.py
@IDE: pycharm
@Description:
    基本回调函数注册表，此时不考虑云端情况
"""
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.data.run.main import get_run, SwanLabRunState
from swanlab.data.run.callback import SwanLabRunCallback
import traceback


class LocalRunCallback(SwanLabRunCallback):

    def __init__(self):
        super(LocalRunCallback, self).__init__()

    @staticmethod
    def _traceback_error(tb):
        """
        获取traceback信息
        """
        trace_list = traceback.format_tb(tb)
        html = ""
        for line in trace_list:
            html += line + "\n"
        return html

    @staticmethod
    def _error_print(tp):
        """
        错误打印
        """
        # 如果是KeyboardInterrupt异常
        if tp == KeyboardInterrupt:
            swanlog.error("KeyboardInterrupt by user")
        else:
            swanlog.error("Error happened while training")

    def before_init_project(self, *args, **kwargs):
        pass

    def _except_handler(self, tp, val, tb):
        """
        异常处理
        """
        self._error_print(tp)
        # 结束运行
        get_run().finish(SwanLabRunState.CRASHED, error=self._traceback_error(tb))
        if tp != KeyboardInterrupt:
            raise tp(val)

    def _clean_handler(self):
        run = get_run()
        if run is None:
            return swanlog.debug("SwanLab Runtime has been cleaned manually.")
        self._train_finish_print()
        # 如果正在运行
        run.finish() if run.is_running else swanlog.debug("Duplicate finish, ignore it.")

    def on_train_begin(self, *args, **kwargs):
        """
        训练开始，注册系统回调
        """
        # 注入系统回调
        self._register_sys_callback()
        # 打印信息
        self._train_begin_print()
        swanlog.info("Experiment_name: " + FONT.yellow(self.settings.exp_name))
        self._watch_tip_print()

    def on_train_end(self, error: str = None):
        """
        训练结束，取消系统回调
        此函数被`run.finish`调用
        """
        # 打印信息
        self._watch_tip_print()
        # 取消注册系统回调
        self._unregister_sys_callback()

    def on_metric_create(self, *args, **kwargs):
        pass

    def on_column_create(self, *args, **kwargs):
        pass
