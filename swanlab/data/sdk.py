#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 18:00:04
@File: swanlab/data/sdk.py
@IDE: vscode
@Description:
    在此处封装swanlab在日志记录模式下的各种接口
"""
import atexit, sys, traceback, os
from datetime import datetime
from .run import SwanLabRun, register
from typing import Optional
from ..log import swanlog
from .modules import DataType
from typing import Dict


run: Optional["SwanLabRun"] = None
inited: bool = False


def init(
    experiment_name: str = None,
    description: str = None,
    config: dict = None,
    log_level: str = None,
    suffix: str = "timestamp",
) -> SwanLabRun:
    """
    Start a new run to track and log.
    Once you have called this function, you can use 'swanlab.log' to log data to the current run.
    Meanwhile, you can use 'swanlab.finish' to finish the current run and close the current experiment.
    After calling this function, SwanLab will begin to record the console output of the current process, and register a callback function to the exit function.

    Parameters
    ----------
    experiment_name : str, optional
        The experiment name you currently have open. If this parameter is not provided, SwanLab will generate one for you by default.
    description : str, optional
        The experiment description you currently have open, used for a more detailed introduction or labeling of the current experiment.
        If you do not provide this parameter, you can modify it later in the web interface.
    config : dict, optional
        Some experiment parameter configurations that can be displayed on the web interface, such as learning rate, batch size, etc.
    log_level : str, optional
        The log level of the current experiment, the default is 'info', you can choose from 'debug', 'info', 'warning', 'error', 'critical'.
    suffix : str, optional
        The suffix of the experiment name, used to distinguish experiments with the same name, the format is yyyy-mm-dd_HH-MM-SS.
        If this parameter is not provided, no suffix will be added.
        At present, only 'timestamp' or None is allowed, and other values will be ignored as 'timestamp'.
    """
    global run, inited
    if inited:
        swanlog.warning("You have already initialized a run, the init function will be ignored")
    else:
        run = register(
            experiment_name=experiment_name,
            description=description,
            config=config,
            log_level=log_level,
            suffix=suffix,
        )
        # 注册异常处理函数
        sys.excepthook = __except_handler
        # 注册清理函数
        atexit.register(__clean_handler)
        swanlog.debug("SwanLab Runtime has initialized")
        swanlog.debug("Swanlab will take over all the print information of the terminal from now on")
        swanlog.info("Run data will be saved locally in " + run.settings.exp_dir)
        swanlog.info("Experiment_name: " + run.settings.exp_name)
        swanlog.info("Run `swanlab watch` to view SwanLab Experiment Dashboard")
        inited = True
    return run


def log(data: Dict[str, DataType], step: int = None):
    """
    Log a row of data to the current run.

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
    if not inited:
        raise RuntimeError("You must call swanlab.data.init() before using log()")
    if inited and run is None:
        return swanlog.error("After calling finish(), you can no longer log data to the current experiment")
    return run.log(data, step)


def finish():
    """
    Finish the current run and close the current experiment
    Normally, swanlab will run this function automatically,
    but you can also execute it manually and mark the experiment as 'completed'.
    Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
    """
    global run
    if not inited:
        raise RuntimeError("You must call swanlab.data.init() before using finish()")
    if run is None:
        return swanlog.error("After calling finish(), you can no longer close the current experiment")
    run.success()
    swanlog.setSuccess()
    swanlog.reset_console()
    run = None


# 定义清理函数
def __clean_handler():
    if run is None:
        return swanlog.debug("SwanLab Runtime has been cleaned manually.")
    if not swanlog.isError:
        swanlog.info(
            "The current experiment {} has been completed, SwanLab will close it automatically".format(
                run.settings.exp_name
            )
        )
        run.success()
        swanlog.setSuccess()
        swanlog.reset_console()


# 定义异常处理函数
def __except_handler(tp, val, tb):
    if run is None:
        return swanlog.warning("SwanLab Runtime has been cleaned manually, the exception will be ignored")
    swanlog.error("Error happended while training, SwanLab will throw it")
    # 标记实验失败
    run.fail()
    swanlog.setError()
    # 记录异常信息
    # 追踪信息
    traceList = traceback.format_tb(tb)
    html = repr(tp) + "\n"
    html += repr(val) + "\n"
    for line in traceList:
        html += line + "\n"

    if os.path.exists(run.settings.error_path):
        swanlog.warning("Error log file already exists, append error log to it")
    # 写入日志文件
    with open(run.settings.error_path, "a") as fError:
        print(datetime.now(), file=fError)
        print(html, file=fError)
    # 重置控制台记录器
    swanlog.reset_console()
    raise tp(val)
