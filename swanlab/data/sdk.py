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
from .run import SwanLabRun, SwanLabConfig, register
from typing import Optional, Union
from ..log import swanlog
from .modules import DataType
from typing import Dict
from ..env import init_env, ROOT
from .utils.file import check_dir_and_create, formate_abs_path
from ..db import Project, connect
from ..utils import version_limit


run: Optional["SwanLabRun"] = None
"""Global runtime instance. After the user calls finish(), run will be set to None."""

inited: bool = False
"""Indicates whether init() has been called in the current process."""

config: Optional["SwanLabConfig"] = SwanLabConfig(None)
"""
Allows users to record experiment configurations through swanlab.config.
Before calling the init() function, config cannot be read or written, even if it is a SwanLabConfig object.
After calling the init() function, swanlab.config is equivalent to run.config.
Configuration information synchronization is achieved through class variables.
When the run object is initialized, it will operate on the SwanLabConfig object to write the configuration.
"""


def init(
    experiment_name: str = None,
    description: str = None,
    config: dict = None,
    logdir: str = None,
    suffix: str = "timestamp",
    log_level: str = None,
    logggings: bool = False,
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
    dir : str, optional
        The directory where the log file is stored, the default is current working directory.
        You can also specify a directory to store the log file, whether using an absolute path or a relative path, but you must ensure that the directory exists.
    suffix : str, optional
        The suffix of the experiment name, used to distinguish experiments with the same name, the format is yyyy-mm-dd_HH-MM-SS.
        If this parameter is not provided, no suffix will be added.
        At present, only 'timestamp' or None is allowed, and other values will be ignored as 'timestamp'.
    """
    global run, inited

    if inited:
        swanlog.warning("You have already initialized a run, the init function will be ignored")
        return run
    # 如果传入了logdir，则将logdir设置为环境变量，代表日志文件存放的路径
    if logdir is not None:
        try:
            logdir = check_dir_and_create(logdir)
        except ValueError:
            raise ValueError("logdir must be a str.")
        except IOError:
            raise IOError("logdir must be a path and have Write permission.")
        os.environ[ROOT] = logdir
    # 如果没有传入logdir，则使用默认的logdir, 即当前工作目录，但是需要保证目录存在
    else:
        logdir = os.path.abspath("swanlog")
        try:
            os.makedirs(logdir, exist_ok=True)
            if not os.access(logdir, os.W_OK):
                raise IOError
        except:
            raise IOError("logdir must have Write permission.")
    # 检查logdir内文件的版本，如果<=0.1.4则报错
    version_limit(logdir, mode="init")
    # 初始化环境变量
    init_env()
    # 连接数据库，要求路径必须存在，但是如果数据库文件不存在，会自动创建
    connect(autocreate=True)

    # 初始化项目数据库
    Project.init(os.path.basename(os.getcwd()))
    # 注册实验
    run = register(
        experiment_name=experiment_name,
        description=description,
        config=config,
        log_level=log_level,
        suffix=suffix,
        loggings=logggings,
    )
    # 注册异常处理函数
    sys.excepthook = __except_handler
    # 注册清理函数
    atexit.register(__clean_handler)
    swanlog.debug("SwanLab Runtime has initialized")
    swanlog.debug("Swanlab will take over all the print information of the terminal from now on")
    swanlog.info("Run data will be saved locally in " + formate_abs_path(run.settings.run_dir))
    swanlog.info("Experiment_name: " + run.settings.exp_name)
    swanlog.info("Run `swanlab watch` to view SwanLab Experiment Dashboard")
    inited = True
    return run


def log(data: Dict[str, DataType], step: int = None, logger: Union[bool, dict] = None):
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
    logger : bool or dict, optional
        Whether to print the data to the console, the default is None.
        If you pass a bool, you can specify whether to print the data to the console.
        If you pass a dict, you can specify whether to print the data to the console, the prefix and suffix of the print data, whether to print timestamp.
        Examples1: swanlab.log({"loss": loss}, logger=True)
        Examples2: swanlab.log({"loss": loss}, logger={"open": True, "prefix": "[0/200] ", "subfix": None, "time":False})
    """
    if not inited:
        raise RuntimeError("You must call swanlab.data.init() before using log()")
    if inited and run is None:
        return swanlog.error("After calling finish(), you can no longer log data to the current experiment")

    # logger_dict的参数为open, prefix, subfix
    # open: 是否开启打印, 如果为None，则根据init的loggings结果走; 如果为True，则开启打印;如果为False，则关闭打印。log->logger的优先级高于init->loggings。
    # prefix: 打印的前缀, 必须为str, float or init。
    # subfix: 打印的后缀, 必须为str, float or init。
    logger_dict = {
        "open": None,
        "prefix": "",
        "subfix": "",
        "time": True,
    }

    # 如果传入的是布尔值且是True，则默认开启打印
    if isinstance(logger, bool):
        logger_dict["open"] = logger
    # 如果传入的是字典，默认开启打印，并将字典中的值赋值给logger_dict
    elif isinstance(logger, dict):
        logger_dict["open"] = True
        logger_dict.update(logger)

        # 字典内参数类型检查
        if not isinstance(logger_dict["open"], bool):
            raise ValueError("logger's open must be a bool")
        if not isinstance(logger_dict["prefix"], (str, float, int)):
            raise ValueError("logger's prefix must be a str, float or init")
        if not isinstance(logger_dict["subfix"], (str, float, int)):
            raise ValueError("logger's subfix must be a str, float or init")
        if not isinstance(logger_dict["time"], bool):
            raise ValueError("logger's time must be a bool")

        # 如果传入的参数超过了open, prefix, subfix，则警告
        if len(logger_dict) > 4:
            swanlog.warning(
                "logger's valid key only has 'open', 'prefix' , 'subfix' and 'time', other parameters will not take effect"
            )
    # 如果传入的是None，则默认关闭打印
    elif logger is None:
        logger_dict["open"] = None
    else:
        raise ValueError("loggings must be a bool or a dict")

    swanlog.set_temporary_logging(logger_dict["open"])
    l = run.log(data, step, logger_dict)
    swanlog.reset_temporary_logging()
    return l


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


def __clean_handler():
    """定义清理函数"""
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
    """定义异常处理函数"""
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
