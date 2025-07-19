#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:58:58
@File: swanlab/data/run/main.py
@IDE: vscode
@Description:
    在此处定义SwanLabRun类并导出
"""
import os
from typing import Any, Dict, Optional, List, Tuple

from swanlab.data.modules import DataWrapper, FloatConvertible, Line, Echarts, PyEchartsBase, PyEchartsTable
from swanlab.env import get_mode
from swanlab.formatter import check_key_format
from swanlab.log import swanlog
from swanlab.swanlab_settings import reset_settings, get_settings
from swanlab.toolkit import MediaType
from .config import SwanLabConfig
from .exp import SwanLabExp
from .helper import SwanLabRunOperator, RuntimeInfo, SwanLabRunState, MonitorCron
from .metadata import get_requirements, get_conda, HardwareCollector
from .public import SwanLabPublicConfig
from ..store import get_run_store, reset_run_store

MAX_LIST_LENGTH = 108


class SwanLabRun:
    """
    The SwanLabRun class is used for logging during a single experiment.
    There should be only one instance of the SwanLabRun class for each experiment.
    """

    def __init__(
        self,
        metadata: dict = None,
        monitor_funcs: List[HardwareCollector] = None,
        run_config: Any = None,
        operator: SwanLabRunOperator = None,
    ):
        """
        Initializing the SwanLabRun class involves configuring the settings and initiating other logging processes.

        Parameters
        ----------
        run_config : Any, optional
            实验参数配置，可以在web界面中显示，如学习率、batch size等
            不需要做任何限制，但必须是字典类型，可被json序列化，否则会报错
        operator : SwanLabRunOperator, optional
            实验操作员，用于批量处理回调函数的调用，如果不提供此参数(为None)，则会自动生成一个实例
        """
        if self.is_started():
            raise RuntimeError("SwanLabRun has been initialized")
        global run, config
        run_store = get_run_store()
        # ---------------------------------- 初始化类内参数 ----------------------------------
        operator = operator or SwanLabRunOperator()
        # 0. 下面的参数会在实验结束后进行副作用清理
        self.__operator = operator
        self.__state = SwanLabRunState.RUNNING
        self.__monitor_cron: Optional[MonitorCron] = None
        self.__config: Optional[SwanLabConfig] = None
        # 1. 设置常规参数
        self.__mode = get_mode()
        self.__public = SwanLabPublicConfig()
        self.__operator.before_run(None)
        # 2. 初始化配置
        if run_store.config is not None:
            config.update(run_store.config)
        config.update(run_config)
        # FIXME 不要使用setattr来修改私有变量
        setattr(config, "_SwanLabConfig__on_setter", self.__operator.on_runtime_info_update)
        self.__config = config
        # 3. 初始化实验
        self.__run_id = run_store.run_id
        assert self.__run_id is not None, "Run ID must be set before initializing SwanLabRun"
        self.__operator.before_init_experiment(
            run_store.run_id,
            run_store.run_name,
            run_store.description,
            run_store.run_colors,
        )
        run = self
        operator.on_run()
        self.__exp = SwanLabExp(operator=operator)
        # ---------------------------------- 初始化完成 ----------------------------------
        # 执行__save，必须在on_run之后，因为on_run之前部分的信息还没完全初始化
        getattr(config, "_SwanLabConfig__save")()
        # 运行时信息采集
        metadata, requirements, conda = _get_runtime_info(metadata)
        operator.on_runtime_info_update(RuntimeInfo(requirements=requirements, conda=conda, metadata=metadata))
        # 定时采集系统信息
        # 测试时不开启此功能
        # resume时不开启此功能
        if "PYTEST_VERSION" not in os.environ and run_store.resume == 'never':
            if monitor_funcs is not None and len(monitor_funcs) != 0:
                swanlog.debug("Monitor on.")

                # 定义定时任务函数
                def monitor_func():
                    monitor_info_list = [f() for f in monitor_funcs]
                    # 剔除其中为None的数据
                    for monitor_info in monitor_info_list:
                        if monitor_info is None:
                            swanlog.debug("Hardware info is empty. Skip it.")
                            continue
                        for info in monitor_info:
                            key, name, value, cfg = (
                                info['key'],
                                info['name'],
                                info['value'],
                                info['config'],
                            )
                            v = DataWrapper(key, [Line(value)], reference="TIME")
                            self.__exp.add(
                                data=v,
                                key=key,
                                name=name,
                                column_config=cfg,
                                column_class="SYSTEM",
                                section_type="SYSTEM",
                            )

                self.__monitor_cron = MonitorCron(monitor_func)

    def __cleanup(self, error: str = None):
        """
        停止部分功能，内部清理时调用
        """
        # 1. 停止硬件监控
        if self.__monitor_cron is not None:
            self.__monitor_cron.cancel()
        # 2. 更新状态
        self.__state = SwanLabRunState.SUCCESS if error is None else SwanLabRunState.CRASHED
        # 3. 触发回调
        if get_settings().log_proxy_type not in ['stderr', 'all']:
            error = None
        self.__operator.on_stop(error)
        # 4. 更新实验 config
        _config = SwanLabConfig(config)
        self.__config = _config
        config.clean()

    def __str__(self) -> str:
        """此类的字符串表示"""
        return "SwanLabRun(run_id={}, mode={}, state={})".format(
            get_run_store().run_id,
            self.__mode,
            self.__state,
        )

    @property
    def id(self) -> Optional[str]:
        """
        The unique identifier for the current run, which is automatically generated by SwanLab.
        User can also set it manually in the `init` function.
        If mode != 'cloud', this will be None.
        """
        if self.mode == "cloud":
            return self.__run_id
        return None

    @property
    def public(self):
        return self.__public

    @property
    def mode(self) -> str:
        return self.__mode

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
        # 1. 设置实验状态
        # 2. 清理 run 内部副作用，并触发 on_stop 回调
        # 4. 重置运行时配置 run_store、重置用户设置 settings
        # 5. 返回old_run
        # 上述步骤中只有客户端对象 client 不清空，其余全局变量全部清理
        if run is None:
            raise RuntimeError("The run object is None, please call `swanlab.init` first.")
        if state == SwanLabRunState.CRASHED and error is None:
            raise ValueError("When the state is 'CRASHED', the error message cannot be None.")
        error = error if state == SwanLabRunState.CRASHED else None
        # 清理内部副作用，触发 on_stop 回调
        getattr(run, "_SwanLabRun__cleanup")(error)
        # 重置输出代理
        swanlog.reset()
        # 清空配置
        reset_run_store()
        reset_settings()
        # 返回旧的 run 对象
        old_run, run = run, None
        return old_run

    @property
    def config(self) -> SwanLabConfig:
        """
        This property allows you to access the 'config' content passed through `init`,
        and allows you to modify it. The latest configuration after each modification
        will be synchronized to the corresponding path by Swanlab. If you have
        enabled the web service, you will notice the changes after refreshing.
        """
        return self.__config

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
            The value must be a `float`, `float convertible object`, `int`, `swanlab.data.BaseType`, or a nested dict.
            For nested dicts, keys will be joined with dots (e.g., {'a': {'b': 1}} becomes {'a.b': 1}).
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
        self.__operator.on_log(data=data, step=step)

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

        # 展平嵌套字典
        flattened_data = _flatten_dict(data)

        log_return = {}
        # 遍历data，记录data
        for k, v in flattened_data.items():
            _k = k
            k = check_key_format(k, auto_cut=True)
            if k != _k:
                # 超过255字符，截断
                swanlog.warning(f"Key {_k} is too long, cut to 255 characters.")
                if k in flattened_data.keys():
                    raise ValueError(f'tag: Not supported too long Key "{_k}" and auto cut failed')
            # ---------------------------------- 包装数据 ----------------------------------
            # 输入为可转换为float的数据类型
            if isinstance(v, (int, float, FloatConvertible)):
                v = DataWrapper(k, [Line(v)])
            elif isinstance(v, (PyEchartsBase, PyEchartsTable)):
                v = DataWrapper(k, [Echarts(v)])
            # 为Line类型或者MediaType类型
            elif isinstance(v, (Line, MediaType)):
                v = DataWrapper(k, [v])
            # 为List[MediaType]、List[PyEchartsBase]或者List[Line]类型，且长度大于0，且所有元素类型相同
            elif (
                isinstance(v, list)
                and len(v) > 0
                and all([isinstance(i, (Line, MediaType, PyEchartsBase, PyEchartsTable)) for i in v])
                and all([i.__class__ == v[0].__class__ for i in v])
            ):
                if len(v) > MAX_LIST_LENGTH:
                    swanlog.warning(f"List length '{k}' is too long, cut to {MAX_LIST_LENGTH}.")
                    v = v[:MAX_LIST_LENGTH]
                # echarts 类型需要转换
                if isinstance(v[0], (PyEchartsBase, PyEchartsTable)):
                    v = DataWrapper(k, [Echarts(i) for i in v])
                else:
                    v = DataWrapper(k, v)
            else:
                # 其余情况被当作是非法的数据类型，交给Line处理
                v = DataWrapper(k, [Line(v)])
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            metric_info = self.__exp.add(key=k, data=v, step=step)
            log_return[metric_info.column_info.key] = metric_info

        return log_return


def _get_runtime_info(metadata: Optional[dict]) -> Tuple[Optional[dict], Optional[str], Optional[str]]:
    """
    获取当前运行时信息，包括requirements、conda和metadata
    :return: requirements, conda, metadata
    """
    user_settings = get_settings()
    requirements, conda = None, None
    if user_settings.requirements_collect:
        requirements = get_requirements()
    if user_settings.conda_collect:
        conda = get_conda()
    if not user_settings.metadata_collect:
        metadata = None
    return metadata, requirements, conda


def _flatten_dict(d: dict, parent_key='', sep='.') -> dict:
    """Helper method to flatten nested dictionaries with dot notation"""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


run: Optional["SwanLabRun"] = None
"""Global runtime instance. After the user calls finish(), run will be set to None."""
config: Optional["SwanLabConfig"] = SwanLabConfig()  # 全局唯一的config对象，不应该重新赋值
"""Global config instance."""


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


def get_url() -> Optional["str"]:
    """
    Get the url of the current experiment.
    NOTE: return None if the experiment has not been initialized or mode is not 'cloud'.
    """
    global run
    if run is None:
        return None
    return run.public.cloud.experiment_url


def get_project_url() -> Optional["str"]:
    """
    Get the url of the current project.
    NOTE: return None if the experiment has not been initialized or mode is not 'cloud'.
    """
    global run
    if run is None:
        return None
    return run.public.cloud.project_url
