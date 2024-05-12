#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback.py
@IDE: pycharm
@Description:
    回调函数注册抽象类
"""
from typing import Union, Tuple, Optional, Callable, Dict
from swanlab.data.settings import SwanDataSettings
from swanlab.data.modules import DataType
from abc import ABC, abstractmethod
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.env import is_windows
from swanlab.package import get_package_version
import atexit
import sys
import os

NewKeyInfo = Union[None, Tuple[dict, Union[float, DataType], int, int]]
"""
新的key对象、数据类型、步数、行数
为None代表没有添加新的key
"""


class ColumnInfo:
    """
    列信息
    """

    def __init__(
        self,
        key: str,
        namespace: str,
        data_type: str,
        chart_type: str,
        sort: int,
        error: Optional[Dict] = None,
        reference: Optional[str] = None,
        config: Optional[Dict] = None
    ):
        self.key = key
        self.namespace = namespace
        self.data_type = data_type
        self.chart_type = chart_type
        self.error = error
        self.reference = reference if reference is not None else "step"
        self.sort = sort
        self.config = config if config is not None else {}


class SwanLabRunCallback(ABC):
    """
    SwanLabRunCallback抽象类
    """

    def __init__(self):
        self.settings: Optional[SwanDataSettings] = None

    def inject(self, settings: SwanDataSettings):
        """
        为SwanLabRunCallback注入settings等一些依赖，因为实例化可能在SwanLabRun之前发生
        :param settings: SwanDataSettings, 数据配置
        :return:
        """
        self.settings = settings

    @staticmethod
    def formate_abs_path(path: str) -> str:
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
        local_path = FONT.magenta(FONT.bold(self.formate_abs_path(self.settings.run_dir)))
        swanlog.info("Run data will be saved locally in " + local_path)

    def _watch_tip_print(self):
        """
        watch命令提示打印
        """
        swanlog.info(
            "🌟 Run `"
            + FONT.bold("swanlab watch -l {}".format(self.formate_abs_path(self.settings.swanlog_dir)))
            + "` to view SwanLab Experiment Dashboard locally"
        )

    def _train_finish_print(self):
        """
        打印结束信息
        """
        swanlog.info("Experiment {} has completed".format(FONT.yellow(self.settings.exp_name)))

    def _register_sys_callback(self):
        """
        注册系统回调
        """
        sys.excepthook = self._except_handler
        atexit.register(self._clean_handler)

    def _unregister_sys_callback(self):
        """
        注销系统回调
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

    @abstractmethod
    def before_init_project(self, proj_name: str, workspace: str):
        """
        在执行业务逻辑之前调用
        """
        pass

    @abstractmethod
    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        suffix: str,
        setter: Callable[[str, str, str, str], None]
    ):
        """
        在初始化实验之前调用
        """
        pass

    @abstractmethod
    def on_log(self):
        """
        每次执行swanlab.log时调用
        """
        pass

    @abstractmethod
    def on_train_begin(self):
        """
        训练开始时的回调函数
        """
        pass

    @abstractmethod
    def on_train_end(self, error: str = None):
        """
        训练结束时的回调函数
        """
        pass

    @abstractmethod
    def on_metric_create(self, key: str, key_info: NewKeyInfo, static_dir: str):
        """
        指标创建回调函数,新增指标信息时调用
        """
        pass

    @abstractmethod
    def on_column_create(self, column_info: ColumnInfo):
        """
        列创建回调函数,新增列信息时调用
        """
        pass
