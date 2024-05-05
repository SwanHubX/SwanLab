#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback.py
@IDE: pycharm
@Description:
    基本回调函数注册表，此时不考虑云端情况
"""
from typing import Union, Tuple, Optional, Dict
from swanlab.data.modules import DataType
from swanlab.data.settings import SwanDataSettings
from swanlab.log import swanlog
from ..utils.file import formate_abs_path
from swanlab.package import get_package_version
from swanlab.utils.font import FONT

NewKeyInfo = Union[None, Tuple[dict, Union[float, DataType], int, int]]
"""
新的key对象、数据类型、步数、行数
为None代表没有添加新的key
"""


class SwanLabRunCallback:

    def __init__(self):
        self.pool = None
        self.settings: Optional[SwanDataSettings] = None

    def inject(self, settings: SwanDataSettings):
        """
        为SwanLabRunCallback注入settings，因为实例化可能在SwanLabRun之前发生
        :param settings: SwanDataSettings, 数据配置
        :return:
        """
        self.settings = settings

    def _train_begin_print(self):
        """
        训练开始时的打印信息
        """
        swanlog.debug("SwanLab Runtime has initialized")
        swanlog.debug("SwanLab will take over all the print information of the terminal from now on")
        swanlog.info("Tracking run with swanlab version " + get_package_version())
        local_path = FONT.magenta(FONT.bold(formate_abs_path(self.settings.run_dir)))
        swanlog.info("Run data will be saved locally in " + local_path)

    def _command_tip_print(self):
        swanlog.info(
            "🌟 Run `"
            + FONT.bold("swanlab watch -l {}".format(formate_abs_path(self.settings.swanlog_dir)))
            + "` to view SwanLab Experiment Dashboard locally"
        )

    def on_train_begin(self, *args, **kwargs):
        """
        训练开始回调函数
        """
        self._train_begin_print()
        swanlog.info("Experiment_name: " + FONT.yellow(self.settings.exp_name))
        self._command_tip_print()

    def on_metric_create(self, *args, **kwargs):
        pass

    def on_column_create(self, *args, **kwargs):
        pass

    def on_train_end(self, *args, **kwargs):
        pass
