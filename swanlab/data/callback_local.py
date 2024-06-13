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
from swanlab.data.run.callback import SwanLabRunCallback, MetricInfo, RuntimeInfo
from datetime import datetime
import traceback
import json
import os


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
            swanlog.info("KeyboardInterrupt by user")
        else:
            swanlog.info("Error happened while training")

    def __str__(self):
        return "SwanLabLocalRunCallback"

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
        run.finish() if run.running else swanlog.debug("Duplicate finish, ignore it.")

    def on_init(self, proj_name: str, workspace: str, logdir: str = None):
        pass

    def on_run(self):
        swanlog.install(self.settings.console_dir)
        # 注入系统回调
        self._register_sys_callback()
        # 打印信息
        self._train_begin_print()
        swanlog.info("Experiment_name: " + FONT.yellow(self.settings.exp_name))
        self._watch_tip_print()

    def on_runtime_info_update(self, r: RuntimeInfo):
        # 更新运行时信息
        if r.requirements is not None:
            r.requirements.write(self.settings.files_dir)
        if r.metadata is not None:
            r.metadata.write(self.settings.files_dir)
        if r.config is not None:
            r.config.write(self.settings.files_dir)

    def on_metric_create(self, metric_info: MetricInfo):
        # 出现任何错误直接返回
        if metric_info.error:
            return
        # ---------------------------------- 保存指标数据 ----------------------------------

        self.settings.mkdir(os.path.dirname(metric_info.metric_path))
        self.settings.mkdir(os.path.dirname(metric_info.summary_path))
        with open(metric_info.summary_path, "w+", encoding="utf-8") as f:
            json.dump(metric_info.summary, f, ensure_ascii=False)
        with open(metric_info.metric_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(metric_info.metric, ensure_ascii=False) + "\n")

        # ---------------------------------- 保存媒体字节流数据 ----------------------------------
        if metric_info.buffers is None:
            return
        for i, r in enumerate(metric_info.buffers):
            if r is None:
                continue
            # 组合路径
            path = os.path.join(self.settings.media_dir, metric_info.key)
            os.makedirs(path, exist_ok=True)
            # 写入数据
            with open(os.path.join(path, metric_info.metric["data"][i]), "wb") as f:
                f.write(r.getvalue())

    def on_stop(self, error: str = None):
        """
        训练结束，取消系统回调
        此函数被`run.finish`调用
        """
        # 写入错误信息
        if error is not None:
            with open(self.settings.error_path, "a") as fError:
                print(datetime.now(), file=fError)
                print(error, file=fError)
        # 打印信息
        self._watch_tip_print()
        # 取消注册系统回调
        self._unregister_sys_callback()
