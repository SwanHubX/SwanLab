#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/19 16:46
@File: callback.py
@IDE: pycharm
@Description:
    回调函数注册抽象模块
"""
from typing import Optional
from swankit.core import SwanLabSharedSettings
from swankit.callback import SwanKitCallback
from swanlab.log import swanlog
from swankit.log import FONT
from swanlab.env import is_windows
from swanlab.package import get_package_version
import atexit
import sys
import os


class U:
    """
    工具函数类，隔离SwanLabRunCallback回调与其他工具函数
    """

    def __init__(self):
        self.settings: Optional[SwanLabSharedSettings] = None

    @staticmethod
    def fmt_windows_path(path: str) -> str:
        """这主要针对windows环境，输入的绝对路径可能不包含盘符，这里进行补充
        主要是用于打印效果
        如果不是windows环境，直接返回path，相当于没有调用这个函数

        Parameters
        ----------
        path : str
            待转换的路径

        Returns
        -------
        str
            增加了盘符的路径
        """
        if not is_windows():
            return path
        if not os.path.isabs(path):
            return path
        need_add = len(path) < 3 or path[1] != ":"
        # 处理反斜杠, 保证路径的正确性
        path = path.replace("/", "\\")
        if need_add:
            return os.path.join(os.getcwd()[:2], path)
        return path

    def _train_begin_print(self):
        """
        训练开始时的打印信息
        """
        swanlog.debug("SwanLab Runtime has initialized")
        swanlog.debug("SwanLab will take over all the print information of the terminal from now on")
        swanlog.info("Tracking run with swanlab version " + get_package_version())
        local_path = FONT.magenta(FONT.bold(self.fmt_windows_path(self.settings.run_dir)))
        swanlog.info("Run data will be saved locally in " + local_path)

    def _watch_tip_print(self):
        """
        watch命令提示打印
        """
        swanlog.info(
            "🌟 Run `"
            + FONT.bold("swanlab watch -l {}".format(self.fmt_windows_path(self.settings.swanlog_dir)))
            + "` to view SwanLab Experiment Dashboard locally"
        )

    def _train_finish_print(self):
        """
        打印结束信息
        """
        swanlog.info("Experiment {} has completed".format(FONT.yellow(self.settings.exp_name)))


class SwanLabRunCallback(SwanKitCallback, U):
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

    def __str__(self):
        raise NotImplementedError("Please implement this method")
