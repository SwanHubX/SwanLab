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
from typing import Optional
from ..log import swanlog
from .modules import DataType
from typing import Dict
from ..env import init_env, ROOT, is_login, get_user_api_key
from .utils.file import check_dir_and_create, formate_abs_path
from ..db import Project, connect
from ..utils import version_limit, FONT, get_package_version
from ..utils.package import get_host_web
from ..auth import get_exp_token, terminal_login, code_login
from ..error import NotLoginError, ValidationError
import asyncio


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


def login(api_key: str):
    """
    Login to SwanLab Cloud. If you already have logged in, you can use this function to relogin.

    Parameters
    ----------
    api_key : str
        authentication key.
    """
    # 如果已经登录且保存，判断一下当前api_key是否和本地api_key一致，如果一致，直接返回
    # 如果不一致，继续下面的步骤
    if is_login() and api_key == get_user_api_key():
        return
    # 否则进行登录
    code_login(api_key)


def init(
    experiment_name: str = None,
    description: str = None,
    config: dict = None,
    config_file: str = None,
    logdir: str = None,
    suffix: str = "default",
    log_level: str = None,
    # cloud: bool = False,
    # project: str = None,
    # organization: str = None,
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
    config_file: str, optional
        The path to the configuration file, the default is None.
        If you provide this parameter, SwanLab will read the configuration from the file and update 'config' you provide.
        The configuration file must be in the format of json or yaml.
    log_level : str, optional
        The log level of the current experiment, the default is 'info', you can choose from 'debug', 'info', 'warning', 'error', 'critical'.
    logdir : str, optional
        The directory where the log file is stored, the default is current working directory.
        You can also specify a directory to store the log file, whether using an absolute path or a relative path, but you must ensure that the directory exists.
    suffix : str, optional
        The suffix of the experiment name, the default is 'default'.
        If this parameter is 'default', suffix will be '%b%d-%h-%m-%s_<hostname>'(example:'Feb03_14-45-37_windowsX'), which represents the current time.
        example: experiment_name = 'example', suffix = 'default' -> 'example_Feb03_14-45-37_windowsX';
        If this parameter is None, no suffix will be added.
        If this paramter is a string, the suffix will be the string you provided.
        Attention: experiment_name + suffix must be unique, otherwise the experiment will not be created.
    cloud : bool, optional
        Whether to use the cloud mode, the default is False.
        If you use the cloud mode, the log file will be stored in the cloud, which will still be saved locally.
        If you are not using cloud mode, the `project` and `organization` fields are invalid.
    project : str, optional
        The project name of the current experiment, the default is None, which means the current project name is the same as the current working directory.
        If you are using cloud mode, you must provide the project name.
    organization : str, optional
        The organization name of the current experiment, the default is None, which means the log file will be stored in your personal space.
    """
    global run, inited
    # ---------------------------------- 一些变量、格式检查 ----------------------------------
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

    # 如果传入了config_file，则检查config_file是否是一个字符串，以及转换为绝对路径
    if config_file is not None:
        if not isinstance(config_file, str):
            raise ValueError("config_file must be a string")
        if not os.path.isabs(config_file):
            config_file = os.path.abspath(config_file)

    # 检查logdir内文件的版本，如果<=0.1.4则报错
    version_limit(logdir, mode="init")
    # 初始化环境变量
    init_env()

    # ---------------------------------- 用户登录、格式、权限校验 ----------------------------------
    # 1. 如果没有登录，提示登录
    # 2. 如果登录了，发起请求，如果请求失败，重新登录，返回步骤1
    # token = _get_exp_token(cloud=cloud)
    # 连接本地数据库，要求路径必须存在，但是如果数据库文件不存在，会自动创建
    connect(autocreate=True)

    # 初始化项目数据库
    Project.init(os.path.basename(os.getcwd()))
    # 注册实验
    run = register(
        experiment_name=experiment_name,
        description=description,
        config=config,
        config_file=config_file,
        log_level=log_level,
        suffix=suffix,
    )
    # 如果使用云端模式，在此开启其他线程负责同步数据

    # 注册异常处理函数
    sys.excepthook = __except_handler
    # 注册清理函数
    atexit.register(__clean_handler)
    swanlog.debug("SwanLab Runtime has initialized")
    swanlog.debug("SwanLab will take over all the print information of the terminal from now on")
    # 展示相关信息信息
    swanlog.info("Tracking run with swanlab version " + get_package_version())
    swanlog.info("Run data will be saved locally in " + FONT.magenta(FONT.bold(formate_abs_path(run.settings.run_dir))))
    # not cloud and swanlog.info("Experiment_name: " + FONT.yellow(run.settings.exp_name))
    swanlog.info("Experiment_name: " + FONT.yellow(run.settings.exp_name))
    # 云端版本有一些额外的信息展示
    # cloud and swanlog.info("Syncing run " + FONT.yellow(run.settings.exp_name) + " to the cloud")
    swanlog.info(
        "🌟 [Offline Dashboard] Run `"
        + FONT.bold("swanlab watch -l {}".format(formate_abs_path(run.settings.swanlog_dir)))
        + "` to view SwanLab Experiment Dashboard locally"
    )
    # project_url = get_host_web() + "/" + "{project_name}"
    # experiment_url = project_url + "/" + token
    # cloud and swanlog.info("🏠 View project at " + FONT.blue(FONT.underline(project_url)))
    # cloud and swanlog.info("🚀 View run at " + FONT.blue(FONT.underline(experiment_url)))
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

    l = run.log(data, step)
    # swanlog.reset_temporary_logging()
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


def _get_exp_token(cloud: bool = False):
    """获取当前实验的相关信息
    可能包含实验的token、实验的id、用户信息等信息
    无论是否使用cloud模式，此函数都会执行，都会返回token，不使用cloud模式返回None，对于后面代码而言，token如果为None，说明没有登录
    """
    token = None
    if cloud:
        # 登录成功会返回当前实验的token
        while True:
            try:
                token = asyncio.run(get_exp_token())
                break
            except NotLoginError:
                # 如果没有登录，提示登录
                terminal_login()
    return token


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
