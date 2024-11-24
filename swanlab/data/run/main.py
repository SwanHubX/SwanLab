#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:58:58
@File: swanlab/data/run/main.py
@IDE: vscode
@Description:
    在此处定义SwanLabRun类并导出
"""
import random
from datetime import datetime
from enum import Enum
from typing import Callable, Optional, Dict

from swankit.core import SwanLabSharedSettings

from swanlab.data.modules import MediaType, DataWrapper, FloatConvertible, Line
from swanlab.env import get_mode, get_swanlog_dir
from swanlab.log import swanlog
from swanlab.package import get_package_version
from . import namer as N
from .config import SwanLabConfig
from .exp import SwanLabExp
from .helper import SwanLabRunOperator, RuntimeInfo
from .metadata import get_requirements, get_metadata
from .public import SwanLabPublicConfig
from ..formater import check_key_format, check_exp_name_format, check_desc_format

MAX_LIST_LENGTH = 108


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
        run_config=None,
        log_level: str = None,
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
        run_config : Any, optional
            实验参数配置，可以在web界面中显示，如学习率、batch size等
            不需要做任何限制，但必须是字典类型，可被json序列化，否则会报错
        log_level : str, optional
            当前实验的日志等级，默认为 'info'，可以从 'debug' 、'info'、'warning'、'error'、'critical' 中选择
            不区分大小写，如果不提供此参数(为None)，则默认为 'info'
            如果提供的日志等级不在上述范围内，默认改为info
        exp_num : int, optional
            历史实验总数，用于云端颜色与本地颜色的对应
        operator : SwanLabRunOperator, optional
            实验操作员，用于批量处理回调函数的调用，如果不提供此参数(为None)，则会自动生成一个实例
        """
        if self.is_started():
            raise RuntimeError("SwanLabRun has been initialized")

        # ---------------------------------- 初始化类内参数 ----------------------------------
        self.__project_name = project_name
        # 生成一个唯一的id，随机生成一个8位的16进制字符串，小写
        _id = hex(random.randint(0, 2**32 - 1))[2:].zfill(8)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.__run_id = "run-{}-{}".format(timestamp, _id)
        # 操作员初始化
        self.__operator = SwanLabRunOperator() if operator is None else operator
        self.__settings = SwanLabSharedSettings(
            logdir=get_swanlog_dir(),
            run_id=self.__run_id,
            should_save=not self.__operator.disabled,
            version=get_package_version(),
        )
        self.__public = SwanLabPublicConfig(self.__project_name, self.__settings)
        self.__operator.before_run(self.__settings)
        # ---------------------------------- 初始化日志记录器 ----------------------------------
        swanlog.level = self.__check_log_level(log_level)
        # ---------------------------------- 初始化配置 ----------------------------------
        # 如果config是以下几个类别之一，则抛出异常
        if isinstance(run_config, (int, float, str, bool, list, tuple, set)):
            raise TypeError(
                f"config: {run_config} (type: {type(run_config)}) is not a json serialized dict "
                f"(Support type is dict, MutableMapping, omegaconf.DictConfig, Argparse.Namespace), please check it"
            )
        global config
        config.update(run_config)
        setattr(config, "_SwanLabConfig__on_setter", self.__operator.on_runtime_info_update)
        self.__config = config
        # ---------------------------------- 注册实验 ----------------------------------
        self.__exp: SwanLabExp = self.__register_exp(experiment_name, description, num=exp_num)
        # 实验状态标记，如果status不为0，则无法再次调用log方法
        self.__state = SwanLabRunState.RUNNING

        # 动态定义一个方法，用于修改实验状态
        def _(state: SwanLabRunState):
            self.__state = state

        global _change_run_state, run
        _change_run_state = _
        run = self

        # ---------------------------------- 初始化完成 ----------------------------------
        self.__operator.on_run()
        # 执行__save，必须在on_run之后，因为on_run之前部分的信息还没完全初始化
        getattr(config, "_SwanLabConfig__save")()
        # 系统信息采集
        self.__operator.on_runtime_info_update(
            RuntimeInfo(
                requirements=get_requirements(),
                metadata=get_metadata(self.__settings.log_dir),
            )
        )

    @property
    def public(self):
        return self.__public

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

    @staticmethod
    def is_started() -> bool:
        """
        If the experiment has been initialized, return True, otherwise return False.
        """
        return get_run() is not None

    @property
    def crashed(self) -> bool:
        """
        If the experiment is marked as 'CRASHED', return True, otherwise return False.
        """
        return self.__state == SwanLabRunState.CRASHED

    @property
    def success(self) -> bool:
        """
        If the experiment is marked as 'SUCCESS', return True, otherwise return False.
        """
        return self.__state == SwanLabRunState.SUCCESS

    @property
    def running(self) -> bool:
        """
        If the experiment is marked as 'RUNNING', return True, otherwise return False.
        """
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
        global run, config
        # 分为几步
        # 1. 设置数据库实验状态为对应状态
        # 2. 判断是否为云端同步，如果是则开始关闭线程池和同步状态
        # 3. 清空run和config对象，run改为局部变量_run，新建一个config对象，原本的config对象内容转移到新的config对象，全局config被清空
        # 4. 返回_run
        if run is None:
            raise RuntimeError("The run object is None, please call `swanlab.init` first.")
        if state == SwanLabRunState.CRASHED and error is None:
            raise ValueError("When the state is 'CRASHED', the error message cannot be None.")
        _set_run_state(state)
        error = error if state == SwanLabRunState.CRASHED else None
        # 退出回调
        getattr(run, "_SwanLabRun__operator").on_stop(error)
        try:
            swanlog.uninstall()
        except RuntimeError:
            # disabled 模式下没有install，所以会报错
            pass

        # ---------------------------------- 清空config和run ----------------------------------
        _config = SwanLabConfig(config)
        setattr(run, "_SwanLabRun__config", _config)
        config.clean()
        _run, run = run, None

        return _run

    @property
    def config(self) -> SwanLabConfig:
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

        Raises
        ----------
        ValueError:
            Unsupported key names.
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
        for k, v in data.items():
            _k = k
            k = check_key_format(k, auto_cut=True)
            if k != _k:
                # 超过255字符，截断
                swanlog.warning(f"Key {_k} is too long, cut to 255 characters.")
                if k in data.keys():
                    raise ValueError(f'tag: Not supported too long Key "{_k}" and auto cut failed')
            # ---------------------------------- 包装数据 ----------------------------------
            # 输入为可转换为float的数据类型
            if isinstance(v, (int, float, FloatConvertible)):
                v = DataWrapper(k, [Line(v)])
            # 为Line类型或者MediaType类型
            elif isinstance(v, (Line, MediaType)):
                v = DataWrapper(k, [v])
            # 为List[MediaType]或者List[Line]类型，且长度大于0，且所有元素类型相同
            elif (
                isinstance(v, list)
                and len(v) > 0
                and all([isinstance(i, (Line, MediaType)) for i in v])
                and all([i.__class__ == v[0].__class__ for i in v])
            ):
                if len(v) > MAX_LIST_LENGTH:
                    swanlog.warning(f"List length '{k}' is too long, cut to {MAX_LIST_LENGTH}.")
                    v = v[:MAX_LIST_LENGTH]
                v = DataWrapper(k, v)
            else:
                # 其余情况被当作是非法的数据类型，交给Line处理
                v = DataWrapper(k, [Line(v)])
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            metric_info = self.__exp.add(key=k, data=v, step=step)
            self.__operator.on_metric_create(metric_info)
            log_return[metric_info.key] = metric_info

        return log_return

    def __str__(self) -> str:
        """此类的字符串表示"""
        return self.__run_id

    def __register_exp(
        self,
        experiment_name: str = None,
        description: str = None,
        num: int = None,
    ) -> SwanLabExp:
        """
        注册实验，将实验配置写入数据库中，完成实验配置的初始化
        """
        if experiment_name:
            e = check_exp_name_format(experiment_name)
            if experiment_name != e:
                swanlog.warning("The experiment name has been truncated automatically.")
                experiment_name = e
        if description:
            d = check_desc_format(description)
            if description != d:
                swanlog.warning("The description has been truncated automatically.")
                description = d
        experiment_name = N.generate_name(num) if experiment_name is None else experiment_name
        description = "" if description is None else description
        colors = N.generate_colors(num)
        self.__operator.before_init_experiment(self.__run_id, experiment_name, description, num, colors)
        self.__settings.exp_name = experiment_name
        self.__settings.exp_colors = colors
        self.__settings.description = description
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

# 全局唯一的config对象，不应该重新赋值
config: Optional["SwanLabConfig"] = SwanLabConfig()
"""
Global config instance.
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


def get_config() -> Optional["SwanLabConfig"]:
    """
    Get the current config object.
    """
    global config
    return config
