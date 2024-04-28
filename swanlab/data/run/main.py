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
from ..system import get_system_info, get_requirements
from swanlab.log import swanlog
from swanlab.utils.file import (
    check_exp_name_format,
    check_desc_format,
)
from swanlab.db import Experiment, ExistedError, NotExistedError
from swanlab.data.modules import BaseType
from swanlab.data.config import SwanLabConfig
from swanlab.cloud import ThreadPool
from swanlab.utils import FONT, create_time
from swanlab.api import get_http
from swanlab.api.upload import upload_logs
import traceback
import os
import sys
import time
import random
import atexit
import ujson
from enum import Enum
from .exp import SwanLabExp
from datetime import datetime
from typing import Tuple, Callable, Optional, Dict


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
        experiment_name: str = None,
        description: str = None,
        config: dict = None,
        log_level: str = None,
        suffix: str = None,
        exp_num: int = None,
        pool: ThreadPool = None,
        callbacks: Dict[str, Callable] = None
    ):
        """
        Initializing the SwanLabRun class involves configuring the settings and initiating other logging processes.

        Parameters
        ----------
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
        pool : ThreadPool
            线程池对象，用于云端同步，传入此对象，run.cloud将被设置为True
        callbacks : Dict[str, Callable]
            回调函数字典，key为回调函数名，value为回调函数
        """
        global run
        if run is not None:
            raise RuntimeError("SwanLabRun has been initialized")
        # ---------------------------------- 初始化类内参数 ----------------------------------
        # 生成一个唯一的id，随机生成一个8位的16进制字符串，小写
        _id = hex(random.randint(0, 2 ** 32 - 1))[2:].zfill(8)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.__run_id = "run-{}-{}".format(timestamp, _id)
        # 初始化配置
        self.__settings = SwanDataSettings(run_id=self.__run_id)
        # ---------------------------------- 初始化日志记录器 ----------------------------------
        # output、console_dir等内容不依赖于实验名称的设置
        swanlog.install(self.__settings.console_dir, self.__check_log_level(log_level))
        # ---------------------------------- 初始化配置 ----------------------------------
        # 给外部1个config
        self.__config = SwanLabConfig(config, self.__settings)
        # ---------------------------------- 注册实验 ----------------------------------
        # 校验描述格式
        description = self.__check_description(description)
        self.__exp: SwanLabExp = self.__register_exp(experiment_name, description, suffix, num=exp_num)
        # 实验状态标记，如果status不为0，则无法再次调用log方法
        self.__state = SwanLabRunState.RUNNING
        self.__pool = pool
        # 设置回调函数
        callbacks is not None and [setattr(self.__exp, key, call(self.pool)) for key, call in callbacks.items()]

        # 动态定义一个方法，用于修改实验状态
        def _(state: SwanLabRunState):
            self.__state = state

        global _change_run_state
        _change_run_state = _
        run = self

    @property
    def cloud(self) -> bool:
        return self.__pool is not None

    @property
    def pool(self) -> ThreadPool:
        return self.__pool

    @property
    def state(self) -> SwanLabRunState:
        return self.__state

    @staticmethod
    def get_state() -> SwanLabRunState:
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
        # 写入错误信息
        if state == SwanLabRunState.CRASHED:
            with open(run.settings.error_path, "a") as fError:
                print(datetime.now(), file=fError)
                print(error, file=fError)
        else:
            error = None
        # 触发云端退出
        if run.cloud and not exiting_in_cloud:
            _before_exit_in_cloud(state, error=error)
        # 重置控制台记录器
        swanlog.uninstall()
        _run, run = run, None
        # 删除回调
        atexit.unregister(clean_handler)
        sys.excepthook = sys.__excepthook__
        return _run

    def set_metric_callback(self, callback: Callable):
        """
        设置实验指标回调函数，当新指标出现时，会调用此函数，只能设置一次，否则出现ValueError
        :param callback: 回调函数

        :raises: ValueError - 如果已经设置过回调函数，再次设置会抛出异常
        """
        self.__exp.metric_callback = callback

    def set_column_callback(self, callback: Callable):
        """
        设置实验列回调函数，当新列出现时，会调用此函数，只能设置一次，否则出现ValueError
        :param callback: 回调函数

        :raises: ValueError - 如果已经设置过回调函数，再次设置会抛出异常
        """
        self.__exp.column_callback = callback

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
        # 每一次log的时候检查一下数据库中的实验状态
        # 如果实验状态不为0，说明实验已经结束，不允许再次调用log方法
        # 这意味着每次log都会进行查询，比较消耗性能，后续考虑采用多进程共享内存的方式进行优化
        swanlog.debug(f"Check experiment and state...")
        try:
            exp = Experiment.get(self.__exp.id)
        except NotExistedError:
            raise KeyboardInterrupt("The experiment has been deleted by the user")
        # 此时self.__state == 0，说明是前端主动停止的
        if exp.status != 0:
            raise KeyboardInterrupt("The experiment has been stopped by the user")

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
        # 遍历data，记录data
        for key in data:
            # 遍历字典的key，记录到本地文件中
            d = data[key]
            # 如果d的数据类型是list，且里面的数据全部为Image类型，则需要转换一下
            if isinstance(d, list) and all([isinstance(i, BaseType) for i in d]) and len(d) > 0:
                # 将d作为输入，构造一个与d相同类型的实例
                d = d[0].__class__(d)
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            self.__exp.add(key=key, data=d, step=step)

    def __str__(self) -> str:
        """此类的字符串表示"""
        return self.__run_id

    @staticmethod
    def __get_exp_name(experiment_name: str = None, suffix: str = None) -> Tuple[str, str]:
        """
        预处理实验名称，如果实验名称过长，截断

        Parameters
        ----------
        experiment_name : str
            实验名称
        suffix : str
            实验名称后缀添加方式，可以为None、"default"或自由后缀，前者代表不添加后缀，后者代表添加时间戳后缀

        Returns
        ----------
        experiment_name : str
            校验后的实验名称
        exp_name : str
            最终的实验名称
        """
        # ---------------------------------- 校验实验名称 ----------------------------------
        experiment_name = "exp" if experiment_name is None else experiment_name
        # 校验实验名称
        experiment_name_checked = check_exp_name_format(experiment_name)
        # 如果实验名称太长的tip
        tip = "The experiment name you provided is too long, it has been truncated automatically."

        # 如果前后长度不一样，说明实验名称被截断了，提醒
        if len(experiment_name_checked) != len(experiment_name) and suffix is None:
            swanlog.warning(tip)

        # 如果suffix为None, 则不添加后缀，直接返回
        if suffix is None or suffix is False:
            return experiment_name_checked, experiment_name

        # 如果suffix为True, 则添加默认后缀
        if suffix is True:
            suffix = "default"

        # suffix必须是字符串
        if not isinstance(suffix, str):
            raise TypeError("The suffix must be a string, but got {}".format(type(suffix)))

        # 如果suffix_checked为default，则设置为默认后缀
        if suffix.lower().strip() == "default":
            # 添加默认后缀
            default_suffix = "{}".format(datetime.now().strftime("%b%d_%H-%M-%S"))
            exp_name = "{}_{}".format(experiment_name_checked, default_suffix)
        else:
            exp_name = "{}_{}".format(experiment_name_checked, suffix)

        # 校验实验名称，如果实验名称过长，截断
        experiment_name_checked = check_exp_name_format(exp_name)
        if len(experiment_name_checked) != len(exp_name):
            swanlog.warning(tip)

        return experiment_name_checked, exp_name

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
        # 这个循环的目的是如果创建失败则等零点五秒重新生成后缀重新创建，直到创建成功
        # 但是由于需要考虑suffix为none不生成后缀的情况，所以需要在except中判断一下
        old = experiment_name
        exp = None
        while True:
            experiment_name, exp_name = self.__get_exp_name(old, suffix)
            try:
                # 获得数据库实例
                exp = Experiment.create(name=exp_name, run_id=self.__run_id, description=description, num=num)
                break
            except ExistedError:
                # 如果suffix名为default，说明是自动生成的后缀，需要重新生成后缀
                if isinstance(suffix, str) and suffix.lower().strip() == "default":
                    swanlog.debug(f"Experiment {exp_name} has existed, try another name...")
                    time.sleep(0.5)
                    continue
                # 其他情况下，说明是用户自定义的后缀，需要报错
                else:
                    Experiment.purely_delete(run_id=self.__run_id)
                    raise ExistedError(f"Experiment {exp_name} has existed in local, please try another name.")

        # 实验创建成功，设置实验相关信息
        self.__settings.exp_name = exp_name
        self.settings.exp_colors = (exp.light, exp.dark)
        self.settings.description = description
        # 执行一些记录操作
        self.__record_exp_config()  # 记录实验配置
        return SwanLabExp(self.__settings, exp.id, exp=exp)

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

    @staticmethod
    def __check_description(description: str) -> str:
        """检查实验描述是否合法"""
        if description is None:
            return ""
        desc = check_desc_format(description)
        if desc != description:
            swanlog.warning("The description has been truncated automatically.")
        return desc

    def __record_exp_config(self):
        """创建实验配置目录 files
        - 创建 files 目录
        - 将实验环境写入 files/swanlab-metadata.json 中
        - 将实验依赖写入 files/requirements.txt 中
        """
        requirements_path = self.__settings.requirements_path
        metadata_path = self.__settings.metadata_path
        # 将实验依赖存入 requirements.txt
        with open(requirements_path, "w") as f:
            f.write(get_requirements())
        # 将实验环境(硬件信息、git信息等等)存入 swanlab-metadata.json
        with open(metadata_path, "w") as f:
            ujson.dump(get_system_info(self.__settings), f)


run: Optional["SwanLabRun"] = None
"""Global runtime instance. After the user calls finish(), run will be set to None."""
exiting_in_cloud = False
"""
Indicates whether the program is exiting in the cloud environment.
"""
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
    # 更新数据库中的实验状态
    run.exp.db.update_status(state.value)


def get_run() -> Optional["SwanLabRun"]:
    """
    Get the current run object. If the experiment has not been initialized, return None.
    """
    global run
    return run


def _before_exit_in_cloud(state: SwanLabRunState, error: str = None):
    """
    在云端环境下，退出之前的处理，需要依次执行线程池中的回调

    Parameters
    ----------
    state : SwanLabRunState
        实验状态
    """
    global exiting_in_cloud
    # 如果正在退出或者run对象为None或者不在云端环境下
    if exiting_in_cloud or run is None or not run.cloud:
        return
    # 标志已经退出（需要在下面的逻辑之前标志）
    exiting_in_cloud = True
    sys.excepthook = except_handler

    def _():
        # 关闭线程池，等待上传线程完成
        run.pool.finish()
        # 上传错误日志
        if error is not None:
            msg = [{"message": error, "create_time": create_time(), "epoch": swanlog.epoch + 1}]
            upload_logs(msg, level="ERROR")

    FONT.loading("Waiting for uploading complete", _)
    get_http().update_state(state == SwanLabRunState.SUCCESS)
    exiting_in_cloud = False
    return


def clean_handler():
    """
    程序执行完毕后，清理资源，用于用户没有调用finish的情况
    """
    if run is None:
        return swanlog.debug("SwanLab Runtime has been cleaned manually.")
    # 如果正在运行，且当前没有正在退出云端环境
    if run.is_running and not exiting_in_cloud:
        swanlog.info("Experiment {} has completed".format(FONT.yellow(run.settings.exp_name)))
        run.finish(SwanLabRunState.SUCCESS)
    else:
        swanlog.debug("Duplicate finish, ignore it.")


def except_handler(tp, val, tb):
    """
    定义异常处理函数，用于程序异常退出时的处理
    此函数触发在clean_handler之前
    """
    if run is None:
        return swanlog.debug("SwanLab Runtime has been cleaned manually.")
    if exiting_in_cloud:
        # FIXME not a good way to fix '\n' problem
        print("")
        swanlog.error("Aborted uploading by user")
        sys.exit(1)
    # 如果是KeyboardInterrupt异常
    if tp == KeyboardInterrupt:
        swanlog.error("KeyboardInterrupt by user")
    else:
        swanlog.error("Error happened while training")
    # 追踪信息
    trace_list = traceback.format_tb(tb)
    html = repr(tp) + "\n"
    html += repr(val) + "\n"
    for line in trace_list:
        html += line + "\n"
    if os.path.exists(run.settings.error_path):
        swanlog.warning("Error log file already exists, append error log to it")
    # 标记实验失败
    run.finish(SwanLabRunState.CRASHED, error=html)
    if tp != KeyboardInterrupt:
        raise tp(val)
