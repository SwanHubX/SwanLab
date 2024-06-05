#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback.py
@IDE: pycharm
@Description:
    回调函数注册抽象类
"""
from typing import Union, Optional, Callable, Dict, List
from abc import ABC, abstractmethod
from swanlab.data.settings import SwanDataSettings
from swanlab.data.modules import ChartType, ErrorInfo, MediaBuffer
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.env import is_windows
from swanlab.package import get_package_version
from urllib.parse import quote
import atexit
import sys
import os


class RuntimeInfo:
    """
    运行时信息，包括系统信息，依赖信息等
    """

    def __init__(self, requirements: str = None, system: dict = None, config: dict = None):
        """
        :param requirements: python依赖
        :param system: 系统信息
        :param config: 上传的配置信息
        """
        pass


class ColumnInfo:
    """
    列信息，当创建列时，会生成这个对象
    """

    def __init__(
        self,
        key: str,
        namespace: str,
        chart: ChartType,
        sort: Optional[int] = None,
        error: Optional[ErrorInfo] = None,
        reference: Optional[str] = None,
        config: Optional[Dict] = None,
    ):
        self.key = key
        """
        列的key名称
        """
        self.namespace = namespace
        """
        列的命名空间
        """
        self.chart = chart
        """
        列的图表类型
        """
        self.error = error
        """
        列的类型错误信息
        """
        self.reference = reference
        """
        列的参考对象
        """
        self.sort = sort
        """
        列在namespace中的排序
        """
        self.config = config if config is not None else {}
        """
        列的额外配置信息
        """

    @property
    def got(self):
        """
        传入的错误类型，如果列出错，返回错误类型，如果没出错，`暂时`返回None
        """
        if self.error is None:
            return None
        return self.error.got

    @property
    def expected(self):
        """
        期望的类型，如果列出错，返回期望的类型，如果没出错，`暂时`返回None
        """
        if self.error is None:
            return None
        return self.error.expected


class MetricInfo:
    """
    指标信息，当新的指标被log时，会生成这个对象
    """
    __SUMMARY_NAME = "_summary.json"

    def __init__(
        self,
        key: str,
        column_info: ColumnInfo,
        error: Optional[ErrorInfo],
        metric: Union[Dict, None] = None,
        summary: Union[Dict, None] = None,
        step: int = None,
        epoch: int = None,
        logdir: str = None,
        metric_file_name: str = None,
        media_dir: str = None,
        buffers: List[MediaBuffer] = None,
    ):
        self.__error = error

        self.key = quote(key, safe="")
        """
        指标的key名称，被quote编码
        """
        self.column_info = column_info
        """
        指标对应的列信息
        """
        self.metric = metric
        """
        指标信息，error时为None
        """
        self.summary = summary
        """
        摘要信息，error时为None
        """
        self.step = step
        """
        当前指标的步数，error时为None
        """
        self.epoch = epoch
        """
        当前指标对应本地的行数，error时为None
        """
        self.metric_path = None if self.error else os.path.join(logdir, self.key, metric_file_name)
        """
        指标文件的路径，error时为None
        """
        self.summary_path = None if self.error else os.path.join(logdir, self.key, self.__SUMMARY_NAME)
        """
        摘要文件的路径，error时为None
        """
        self.media_dir = media_dir
        """
        静态文件的根文件夹
        """
        self.buffers = buffers
        """
        需要上传的媒体数据，比特流，error时为None，如果上传为非媒体类型（或Text类型），也为None
        """
        # 写入文件名称，对应上传时的文件名称
        if self.buffers is not None:
            for i, buffer in enumerate(self.buffers):
                buffer.file_name = "{}/{}".format(self.key, metric["data"][i])

    @property
    def error(self) -> bool:
        """
        这条指标信息是否有错误，错误分几种：
            1. 列错误，列一开始就出现问题
            2. 重复错误
            3. 指标错误
        """
        return self.error_info is not None or self.column_error_info is not None

    @property
    def column_error_info(self) -> Optional[ErrorInfo]:
        """
        列错误信息
        """
        return self.column_info.error

    @property
    def error_info(self) -> Optional[ErrorInfo]:
        """
        指标错误信息
        """
        return self.__error

    @property
    def duplicated_error(self) -> bool:
        """
        是否是重复的指标
        """
        return self.__error and self.__error.duplicated

    @property
    def data(self) -> Union[Dict, None]:
        """
        指标数据的data字段
        """
        if self.error:
            return None
        return self.metric["data"]


class ErrorInfo:
    """
    错误信息，当发生错误时，会生成这个对象
    """

    def __init__(
        self,
        desc: str,
    ):
        self.desc = desc
        """
        错误描述
        """

    def __str__(self):
        return f"SwanLabError: {self.desc}"


class U:
    """
    工具函数类，隔离SwanLabRunCallback回调与其他工具函数
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
    def formate_windows_path(path: str) -> str:
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
        local_path = FONT.magenta(FONT.bold(self.formate_windows_path(self.settings.run_dir)))
        swanlog.info("Run data will be saved locally in " + local_path)

    def _watch_tip_print(self):
        """
        watch命令提示打印
        """
        swanlog.info(
            "🌟 Run `"
            + FONT.bold("swanlab watch -l {}".format(self.formate_windows_path(self.settings.swanlog_dir)))
            + "` to view SwanLab Experiment Dashboard locally"
        )

    def _train_finish_print(self):
        """
        打印结束信息
        """
        swanlog.info("Experiment {} has completed".format(FONT.yellow(self.settings.exp_name)))


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

    def on_run_error_from_operator(self, e: ErrorInfo):
        """
        SwanLabRun初始化错误时被操作员调用
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
