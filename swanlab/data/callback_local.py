#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback_local.py
@IDE: pycharm
@Description:
    基本回调函数注册表，此时不考虑云端情况
"""
from typing import Callable
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.data.run.main import get_run, SwanLabRunState
from swanlab.data.run.callback import SwanLabRunCallback, MetricInfo
from swanlab.data.system import get_system_info, get_requirements
from swanlab.env import ROOT
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
            swanlog.error("KeyboardInterrupt by user")
        else:
            swanlog.error("Error happened while training")

    @staticmethod
    def _init_logdir(logdir: str = None) -> str:
        """
        根据传入的logdir，初始化日志文件夹
        FIXME shit code
        """
        # 如果传入了logdir，则将logdir设置为环境变量，代表日志文件存放的路径
        if logdir is not None:
            try:
                if not isinstance(logdir, str):
                    raise ValueError("path must be a string")
                if not os.path.isabs(logdir):
                    logdir = os.path.abspath(logdir)
                # 如果创建失败，也是抛出IOError
                try:
                    os.makedirs(logdir, exist_ok=True)
                except Exception as e:
                    raise IOError(f"create path: {logdir} failed, error: {e}")
                if not os.access(logdir, os.W_OK):
                    raise IOError(f"no write permission for path: {logdir}")
            except ValueError:
                raise ValueError("logdir must be a str.")
            except IOError:
                raise IOError("logdir must be a path and have Write permission.")
            os.environ[ROOT] = logdir
        # 如果没有传入logdir，则使用默认的logdir, 即当前工作目录下的swanlog文件夹，但是需要保证目录存在
        else:
            logdir = os.environ.get(ROOT) or os.path.join(os.getcwd(), "swanlog")
            logdir = os.path.abspath(logdir)
            try:
                os.makedirs(logdir, exist_ok=True)
                if not os.access(logdir, os.W_OK):
                    raise IOError
            except IOError:
                raise IOError("logdir must have Write permission.")
        # 如果logdir是空的，创建.gitignore文件，写入*
        if not os.listdir(logdir):
            with open(os.path.join(logdir, ".gitignore"), "w") as f:
                f.write("*")
        return logdir

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
        run.finish() if run.is_running else swanlog.debug("Duplicate finish, ignore it.")

    def on_init(self, proj_name: str, workspace: str, logdir: str = None):
        self._init_logdir(logdir)

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        suffix: str,
        setter: Callable[[str, str, str, str], None]
    ):
        requirements_path = self.settings.requirements_path
        metadata_path = self.settings.metadata_path
        # 将实验依赖存入 requirements.txt
        with open(requirements_path, "w") as f:
            f.write(get_requirements())
        # 将实验环境(硬件信息、git信息等等)存入 swanlab-metadata.json
        with open(metadata_path, "w") as f:
            json.dump(get_system_info(self.settings), f)

    def on_run(self):
        swanlog.install(self.settings.console_dir)
        # 注入系统回调
        self._register_sys_callback()
        # 打印信息
        self._train_begin_print()
        swanlog.info("Experiment_name: " + FONT.yellow(self.settings.exp_name))
        self._watch_tip_print()
        if not os.path.exists(self.settings.log_dir):
            os.mkdir(self.settings.log_dir)

    def on_metric_create(self, metric_info: MetricInfo):
        if metric_info.error:
            return
        self.settings.mkdir(os.path.dirname(metric_info.metric_path))
        self.settings.mkdir(os.path.dirname(metric_info.summary_path))
        with open(metric_info.summary_path, "w+") as f:
            json.dump(metric_info.summary, f, ensure_ascii=False)
        with open(metric_info.metric_path, "a") as f:
            f.write(json.dumps(metric_info.metric, ensure_ascii=False) + "\n")

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
