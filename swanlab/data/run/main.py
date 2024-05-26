#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:58:58
@File: swanlab/data/run/main.py
@IDE: vscode
@Description:
    在此处定义SwanLabRun类并导出
"""
from ..settings import SwanDataSettings
from swanlab.log import swanlog
from swanlab.data.modules import BaseType
from swanlab.data.config import SwanLabConfig
import random
from enum import Enum
from .exp import SwanLabExp
from datetime import datetime
from typing import Callable, Optional, Dict
from .operator import SwanLabRunOperator
from swanlab.env import get_mode, SwanLabMode


class SwanLabRunState(Enum):
    """SwanLabRunState is an enumeration class that represents the state of the experiment.
    We Recommend that you use this enumeration class to represent the state of the experiment.
    """
    NOT_STARTED = -2
    SUCCESS = 1
    CRASHED = -1
    RUNNING = 0


class SwanLabRun:
    """
    The SwanLabRun class is used for logging during a single experiment.
    There should be only one instance of the SwanLabRun class for each experiment.
    """

    def __init__(
        self,
        project_name: str = None,
        experiment_name: str = None,
        description: str = None,
        config: dict = None,
        log_level: str = None,
        suffix: str = None,
        exp_num: int = None,
        operator: SwanLabRunOperator = SwanLabRunOperator(),
    ):
        """
        Initializing the SwanLabRun class involves configuring the settings and initiating other logging processes.

        Parameters
        ----------
        project_name : str, optional
            项目名称，目前单纯做一个记录
        experiment_name : str, optional
            实验名称，实验名称应该唯一，由0-9，a-z，A-Z，" ","_","-","/"组成
            如果不提供此参数(为None)，SwanLab将自动生成一个实验名称
        description : str, optional
            实验描述，用于对当前实验进行更详细的介绍或标注
            如果不提供此参数(为None)，可以在web界面中进行修改,这意味着必须在此改为空字符串""
        config : dict, optional
            实验参数配置，可以在web界面中显示，如学习率、batch size等
            不需要做任何限制，但必须是字典类型，可被json序列化，否则会报错
        log_level : str, optional
            当前实验的日志等级，默认为 'info'，可以从 'debug' 、'info'、'warning'、'error'、'critical' 中选择
            不区分大小写，如果不提供此参数(为None)，则默认为 'info'
            如果提供的日志等级不在上述范围内，默认改为info
        suffix : str, optional
            实验名称后缀，用于区分同名实验，格式为yyyy-mm-dd_HH-MM-SS
            如果不提供此参数(为None)，不会添加后缀
        exp_num : int, optional
            历史实验总数，用于云端颜色与本地颜色的对应
        operator : SwanLabRunOperator, optional
            实验操作员，用于批量处理回调函数的调用，如果不提供此参数(为None)，则会自动生成一个实例
        """
        global run
        if run is not None:
            raise RuntimeError("SwanLabRun has been initialized")
        # ---------------------------------- 初始化类内参数 ----------------------------------
        self.__project_name = project_name
        # 生成一个唯一的id，随机生成一个8位的16进制字符串，小写
        _id = hex(random.randint(0, 2 ** 32 - 1))[2:].zfill(8)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.__run_id = "run-{}-{}".format(timestamp, _id)
        # 操作员初始化
        self.__operator = SwanLabRunOperator() if operator is None else operator
        self.__settings = SwanDataSettings(run_id=self.__run_id, should_save=not self.__operator.disabled)
        self.__operator.inject(self.__settings)
        # ---------------------------------- 初始化日志记录器 ----------------------------------
        # output、console_dir等内容不依赖于实验名称的设置
        swanlog.set_level(self.__check_log_level(log_level))
        # ---------------------------------- 初始化配置 ----------------------------------
        # 给外部1个config
        self.__config = SwanLabConfig(config, self.__settings)
        # ---------------------------------- 注册实验 ----------------------------------
        self.__exp: SwanLabExp = self.__register_exp(experiment_name, description, suffix, num=exp_num)
        # 实验状态标记，如果status不为0，则无法再次调用log方法
        self.__state = SwanLabRunState.RUNNING

        # 动态定义一个方法，用于修改实验状态
        def _(state: SwanLabRunState):
            self.__state = state

        global _change_run_state
        _change_run_state = _
        run = self

        # ---------------------------------- 初始化完成 ----------------------------------
        self.__operator.on_run()

    @property
    def operator(self) -> SwanLabRunOperator:
        return self.__operator

    @property
    def project_name(self) -> str:
        return self.__project_name

    @property
    def mode(self) -> str:
        return get_mode()

    @property
    def state(self) -> SwanLabRunState:
        return self.__state

    @classmethod
    def get_state(cls) -> SwanLabRunState:
        """
        获取当前实验状态
        """
        global run
        return run.state if run is not None else SwanLabRunState.NOT_STARTED

    @property
    def is_crashed(self) -> bool:
        return self.__state == SwanLabRunState.CRASHED

    @property
    def is_success(self) -> bool:
        return self.__state == SwanLabRunState.SUCCESS

    @property
    def is_running(self) -> bool:
        return self.__state == SwanLabRunState.RUNNING

    @staticmethod
    def finish(state: SwanLabRunState = SwanLabRunState.SUCCESS, error=None):
        """
        Finish the current run and close the current experiment
        Normally, swanlab will run this function automatically,
        but you can also execute it manually and mark the experiment as 'completed'.
        Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
        After calling this function, you can re-run `swanlab.init` to start a new experiment.

        :param state: The state of the experiment, it can be 'SUCCESS', 'CRASHED' or 'RUNNING'.
        :param error: The error message when the experiment is marked as 'CRASHED'. If not 'CRASHED', it should be None.
        """
        global run
        # 分为几步
        # 1. 设置数据库实验状态为对应状态
        # 2. 判断是否为云端同步，如果是则开始关闭线程池和同步状态
        # 3. 清空run对象，run改为局部变量_run
        # 4. 返回_run
        if run is None:
            raise RuntimeError("The run object is None, please call `swanlab.init` first.")
        if state == SwanLabRunState.CRASHED and error is None:
            raise ValueError("When the state is 'CRASHED', the error message cannot be None.")
        _set_run_state(state)
        error = error if state == SwanLabRunState.CRASHED else None
        # 退出回调
        run.operator.on_stop(error)
        try:
            swanlog.uninstall()
        except RuntimeError:
            # disabled 模式下没有install，所以会报错
            pass
        _run, run = run, None
        return _run

    @property
    def settings(self) -> SwanDataSettings:
        """
        This property allows you to access the 'settings' content passed through `init`,
        and runtime settings can not be modified.
        """
        return self.__settings

    @property
    def config(self):
        """
        This property allows you to access the 'config' content passed through `init`,
        and allows you to modify it. The latest configuration after each modification
        will be synchronized to the corresponding path by Swanlab. If you have
        enabled the web service, you will notice the changes after refreshing.
        """
        return self.__config

    @property
    def exp(self) -> SwanLabExp:
        """
        Get the current experiment object. This object is used to log data and control the experiment.
        """
        return self.__exp

    def log(self, data: dict, step: int = None):
        """
        Log a row of data to the current run. Unlike `swanlab.log`, this api will be called directly based on the
        SwanRun instance, removing the initialization process. Of course, after you call the success/fail method,
        this method will be banned from calling.

        Parameters
        ----------
        data : Dict[str, DataType]
            Data must be a dict.
            The key must be a string with 0-9, a-z, A-Z, " ", "_", "-", "/".
            The value must be a `float`, `float convertible object`, `int` or `swanlab.data.BaseType`.
        step : int, optional
            The step number of the current data, if not provided, it will be automatically incremented.
            If step is duplicated, the data will be ignored.
        """
        if self.__state != SwanLabRunState.RUNNING:
            raise RuntimeError("After experiment finished, you can no longer log data to the current experiment")
        self.__operator.on_log()

        if not isinstance(data, dict):
            return swanlog.error(
                "log data must be a dict, but got {}, SwanLab will ignore records it.".format(type(data))
            )
        # 检查step的类型
        if step is not None and (not isinstance(step, int) or step < 0):
            swanlog.error(
                "'step' must be an integer not less than zero, but got {}, SwanLab will automatically set step".format(
                    step
                )
            )
            step = None

        log_return = {}
        # 遍历data，记录data
        for key in data:
            # 遍历字典的key，记录到本地文件中
            d = data[key]
            # 如果d的数据类型是list，且里面的数据全部为Image类型，则需要转换一下
            if isinstance(d, list) and all([isinstance(i, BaseType) for i in d]) and len(d) > 0:
                # 将d作为输入，构造一个与d相同类型的实例
                d = d[0].__class__(d)
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            metric_info = self.__exp.add(key=key, data=d, step=step)
            self.__operator.on_metric_create(metric_info)
            log_return[metric_info.key] = metric_info

        return log_return

    def __str__(self) -> str:
        """此类的字符串表示"""
        return self.__run_id

    def __register_exp(
        self,
        experiment_name: str,
        description: str = None,
        suffix: str = None,
        num: int = None,
    ) -> SwanLabExp:
        """
        注册实验，将实验配置写入数据库中，完成实验配置的初始化
        """

        # ---------------------------------- 初始化实验 ----------------------------------

        def setter(exp_name: str, light_color: str, dark_color: str, desc: str):
            """
            设置实验相关信息
            :param exp_name: 实验名称
            :param light_color: 亮色
            :param dark_color: 暗色
            :param desc: 实验描述
            :return:
            """
            # 实验创建成功，设置实验相关信息
            self.__settings.exp_name = exp_name
            self.settings.exp_colors = (light_color, dark_color)
            self.settings.description = desc

        self.__operator.before_init_experiment(self.__run_id, experiment_name, description, num, suffix, setter)
        return SwanLabExp(self.__settings, operator=self.__operator)

    @staticmethod
    def __check_log_level(log_level: str) -> str:
        """检查日志等级是否合法"""
        valid = ["debug", "info", "warning", "error", "critical"]
        if log_level is None:
            return "info"
        elif log_level.lower() in valid:
            return log_level.lower()
        else:
            swanlog.warning(f"The log level you provided is not valid, it has been set to {log_level}.")
            return "info"


run: Optional["SwanLabRun"] = None
"""Global runtime instance. After the user calls finish(), run will be set to None."""
_change_run_state: Optional["Callable"] = None
"""
修改实验状态的函数，用于在实验状态改变时调用
"""


def _set_run_state(state: SwanLabRunState):
    """
    设置实验状态，只能在run对象存在的情况下调用
    """
    if run is None:
        raise RuntimeError("The run object is None, please call swanlab.init first.")

    if state not in SwanLabRunState or state in [SwanLabRunState.NOT_STARTED, SwanLabRunState.RUNNING]:
        swanlog.warning("Invalid state when set, state must be in SwanLabRunState and not be RUNNING or NOT_STARTED")
        swanlog.warning("SwanLab will set state to `CRASHED`")
        state = SwanLabRunState.CRASHED
    # 设置state
    _change_run_state(state)


def get_run() -> Optional["SwanLabRun"]:
    """
    Get the current run object. If the experiment has not been initialized, return None.
    """
    global run
    return run
