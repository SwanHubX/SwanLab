#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 18:00:04
@File: swanlab/data/sdk.py
@IDE: vscode
@Description:
    在此处封装swanlab在日志记录模式下的各种接口
"""
import atexit
import os
import sys
from typing import Optional, Union, Dict
from .modules import DataType
from .run import (
    SwanLabRunState,
    SwanLabRun,
    register,
    get_run,
    except_handler,
    clean_handler,
)
from .config import SwanLabConfig
from .utils.file import check_dir_and_create, formate_abs_path
from ..db import Project, connect, Experiment
from ..env import init_env, ROOT, get_swanlab_folder
from ..log import swanlog
from ..utils import FONT, check_load_json_yaml, check_proj_name_format
from ..utils.key import get_key
from ..utils.judgment import in_jupyter, show_button_html
from swanlab.api import create_http, get_http, code_login, LoginInfo, terminal_login
from swanlab.api.upload.model import ColumnModel
from swanlab.package import version_limit, get_package_version, get_host_api, get_host_web
from swanlab.error import KeyFileError, ApiError
from swanlab.cloud import LogSnifferTask, ThreadPool
from swanlab.cloud import UploadType

_config: Optional["SwanLabConfig"] = SwanLabConfig(None)
"""
Allows users to record experiment configurations through swanlab.config.
Before calling the init() function, config cannot be read or written, even if it is a SwanLabConfig object.
After calling the init() function, swanlab.config is equivalent to run.config.
Configuration information synchronization is achieved through class variables.
When the run object is initialized, it will operate on the SwanLabConfig object to write the configuration.
"""

login_info = None
"""
Record user login information in the cloud environment.
`swanlab.login` will assign or update this variable.
"""


def _check_proj_name(name: str) -> str:
    """检查项目名称是否合法，如果不合法则抛出ValueError异常
    项目名称必须是一个非空字符串，长度不能超过255个字符

    Parameters
    ----------
    name : str
        待检查的项目名称

    Returns
    -------
    str
        返回项目名称

    Raises
    ------
    ValueError
        项目名称不合法
    """
    _name = check_proj_name_format(name)
    if len(name) != len(_name):
        swanlog.warning(f"project name is too long, auto cut to {_name}")
    return _name


def _create_metric_callback(pool: ThreadPool):
    """
    创建指标回调函数
    """

    def _metric_callback(is_scalar: dict, data):
        """
        上传指标回调函数
        :param is_scalar: 是否为标量指标
        :param data: 上传的指标
        """
        pool.queue.put((UploadType.SCALAR_METRIC if is_scalar else UploadType.MEDIA_METRIC, [data]))

    return _metric_callback


def _create_column_callback(pool: ThreadPool):
    """
    创建列回调函数
    """

    def _column_callback(key, data_type: str, error: Optional[Dict] = None):
        """
        上传列信息回调函数
        :param key: 列名
        :param data_type: 列数据类型
        :param error: 错误信息
        """
        pool.queue.put((UploadType.COLUMN, [ColumnModel(key, data_type.upper(), error)]))

    return _column_callback


def _is_inited():
    """检查是否已经初始化"""
    return get_run() is not None


def login(api_key: str = None):
    """
    Login to SwanLab Cloud. If you already have logged in, you can use this function to relogin.
    Every time you call this function, the previous login information will be overwritten.
    [Note that] this function should be called before `init`.

    Parameters
    ----------
    api_key : str
        authentication key, if not provided, the key will be read from the key file.
    """
    if _is_inited():
        raise RuntimeError("You must call swanlab.login() before using init()")
    global login_info
    login_info = code_login(api_key) if api_key else _login_in_init()


def init(
    project: str = None,
    workspace: str = None,
    experiment_name: str = None,
    description: str = None,
    config: Union[dict, str] = None,
    logdir: str = None,
    suffix: Union[str, None, bool] = "default",
    cloud: bool = True,
    load: str = None,
    **kwargs,
) -> SwanLabRun:
    """
    Start a new run to track and log. Once you have called this function, you can use 'swanlab.log' to log data to
    the current run. Meanwhile, you can use 'swanlab.finish' to finish the current run and close the current
    experiment. After calling this function, SwanLab will begin to record the console output of the current process,
    and register a callback function to the exit function.

    Parameters
    ----------
    project : str, optional
        The project name of the current experiment, the default is None,
        which means the current project name is the same as the current working directory.
    workspace : str, optional
        Where the current project is located, it can be an organization or a user (currently only supports yourself).
        The default is None, which means the current entity is the same as the current user.
    experiment_name : str, optional
        The experiment name you currently have open. If this parameter is not provided,
        SwanLab will generate one for you by default.
    description : str, optional
        The experiment description you currently have open,
        used for a more detailed introduction or labeling of the current experiment.
        If you do not provide this parameter, you can modify it later in the web interface.
    config : Union[dict, str], optional
        If you provide as a dict, it will be used as the configuration of the current experiment.
        If you provide as a string, SwanLab will read the configuration from the file.
        And the configuration file must be in the format of `json` or `yaml`.
        Anyway, you can modify the configuration later after this function is called.
    logdir : str, optional
        The folder will store all the log information generated during the execution of SwanLab.
        If the parameter is None,
        SwanLab will generate a folder named "swanlog" in the same path as the code execution to store the data.
        If you want to visualize the generated log files,
        simply run the command `swanlab watch` in the same path where the code is executed
        (without entering the "swanlog" folder).
        You can also specify your own folder, but you must ensure that the folder exists and preferably does not contain
        anything other than data generated by Swanlab.
        In this case, if you want to view the logs,
        you must use something like `swanlab watch -l ./your_specified_folder` to specify the folder path.
    suffix : str, optional
        The suffix of the experiment name, the default is 'default'.
        If this parameter is 'default', suffix will be '%b%d-%h-%m-%s'(example:'Feb03_14-45-37'),
        which represents the current time.
        example: experiment_name = 'example', suffix = 'default' -> 'example_Feb03_14-45-37';
        If this parameter is None or False, no suffix will be added.
        If this parameter is a string, the suffix will be the string you provided.
        Attention: experiment_name + suffix must be unique, otherwise the experiment will not be created.
    cloud : bool, optional
        Whether to use the cloud mode, the default is True.
        If you use the cloud mode, the log file will be stored in the cloud, which will still be saved locally.
        If you are not using cloud mode, the `workspace` fields are invalid.
    load : str, optional
        If you pass this parameter,SwanLab will search for the configuration file you specified
        (which must be in JSON or YAML format)
        and automatically fill in some explicit parameters of this function for you
        (excluding parameters in `**kwargs` and the parameters if they are None).
        In terms of priority, if the parameters passed to init are `None`,
        SwanLab will attempt to replace them from the configuration file you provided;
        otherwise, it will use the parameters you passed as the definitive ones.
    """
    # 如果已经初始化过了，直接返回run
    run = get_run()
    if run is not None:
        swanlog.warning("You have already initialized a run, the init function will be ignored")
        return run
    # ---------------------------------- 一些变量、格式检查 ----------------------------------
    # 如果传入了load，则加载load文件，如果load文件不存在，报错
    if load:
        load_data = check_load_json_yaml(load, load)
        # 尝试更改传入的参数为None的情况，如果传入的参数不是None，不做任何操作
        experiment_name = _load_data(load_data, "experiment_name", experiment_name)
        description = _load_data(load_data, "description", description)
        config = _load_data(load_data, "config", config)
        logdir = _load_data(load_data, "logdir", logdir)
        suffix = _load_data(load_data, "suffix", suffix)
        cloud = _load_data(load_data, "cloud", cloud)
        project = _load_data(load_data, "project", project)
        workspace = _load_data(load_data, "workspace", workspace)
    # 默认实验名称为当前目录名
    project = _check_proj_name(project if project else os.path.basename(os.getcwd()))
    # 初始化logdir参数，接下来logdir被设置为绝对路径且当前程序有写权限
    logdir = _init_logdir(logdir)
    # 初始化confi参数
    config = _init_config(config)
    # 检查logdir内文件的版本，如果<=0.1.4则报错
    version_limit(logdir, mode="init")
    # 初始化环境变量
    init_env()
    # 历史实验总数
    exp_num = None
    # ---------------------------------- 用户登录、格式、权限校验 ----------------------------------
    global login_info
    http, pool = None, None
    if login_info is None and cloud:
        # 用户登录
        login_info = _login_in_init()
    if cloud:
        http = create_http(login_info)
        exp_num = http.mount_project(project, workspace).history_exp_count
        # 初始化、挂载线程池
        pool = ThreadPool()
    # 连接本地数据库，要求路径必须存在，但是如果数据库文件不存在，会自动创建
    connect(autocreate=True)
    # 初始化项目数据库
    Project.init(project)
    # ---------------------------------- 实例化实验 ----------------------------------
    # 如果是云端环境，设置回调函数
    callbacks = (
        None if not cloud else {"metric_callback": _create_metric_callback, "column_callback": _create_column_callback}
    )

    # 注册实验
    run = register(
        experiment_name=experiment_name,
        description=description,
        config=config,
        log_level=kwargs.get("log_level", "info"),
        suffix=suffix,
        exp_num=exp_num,
        pool=pool,
        callbacks=callbacks,
    )
    # ---------------------------------- 注册实验，开启线程 ----------------------------------
    if cloud:
        # 注册实验信息
        try:
            get_http().mount_exp(
                exp_name=run.settings.exp_name,
                colors=run.settings.exp_colors,
                description=run.settings.description,
            )
        except ApiError as e:
            if e.resp.status_code == 409:
                FONT.brush("", 50)
                swanlog.error("The experiment name already exists, please change the experiment name")
                Experiment.purely_delete(run_id=run.exp.db.run_id)
                sys.exit(409)
        sniffer = LogSnifferTask(run.settings.files_dir)
        run.pool.create_thread(sniffer.task, name="sniffer", callback=sniffer.callback)

        def _write_call_call(message):
            pool.queue.put((UploadType.LOG, [message]))

        swanlog.set_write_callback(_write_call_call)
    # ---------------------------------- 异常处理、程序清理 ----------------------------------
    sys.excepthook = except_handler
    # 注册清理函数
    atexit.register(clean_handler)
    # ---------------------------------- 终端输出 ----------------------------------
    if not cloud and workspace is not None:
        swanlog.warning("The `workspace` field is invalid in local mode")
    swanlog.debug("SwanLab Runtime has initialized")
    swanlog.debug("SwanLab will take over all the print information of the terminal from now on")
    swanlog.info("Tracking run with swanlab version " + get_package_version())
    swanlog.info("Run data will be saved locally in " + FONT.magenta(FONT.bold(formate_abs_path(run.settings.run_dir))))
    not cloud and swanlog.info("Experiment_name: " + FONT.yellow(run.settings.exp_name))
    # 云端版本有一些额外的信息展示
    cloud and swanlog.info("👋 Hi " + FONT.bold(FONT.default(login_info.username)) + ", welcome to swanlab!")
    cloud and swanlog.info("Syncing run " + FONT.yellow(run.settings.exp_name) + " to the cloud")
    swanlog.info(
        "🌟 Run `"
        + FONT.bold("swanlab watch -l {}".format(formate_abs_path(run.settings.swanlog_dir)))
        + "` to view SwanLab Experiment Dashboard locally"
    )
    if cloud:
        project_url = get_host_web() + f"/@{http.groupname}/{http.projname}"
        experiment_url = project_url + f"/runs/{http.exp_id}"
        swanlog.info("🏠 View project at " + FONT.blue(FONT.underline(project_url)))
        swanlog.info("🚀 View run at " + FONT.blue(FONT.underline(experiment_url)))

        # 在Jupyter Notebook环境下，显示按钮
        if in_jupyter():
            show_button_html(experiment_url)

    return run


def log(data: Dict[str, DataType], step: int = None):
    """
    Log a row of data to the current run.
    We recommend that you log data by SwanLabRun.log() method, but you can also use this function to log data.

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
    if not _is_inited():
        raise RuntimeError("You must call swanlab.init() before using log()")
    run = get_run()
    ll = run.log(data, step)
    return ll


def finish(state: SwanLabRunState = SwanLabRunState.SUCCESS, error=None):
    """
    Finish the current run and close the current experiment
    Normally, swanlab will run this function automatically,
    but you can also execute it manually and mark the experiment as 'completed'.
    Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
    If you mark the experiment as 'CRASHED' manually, `error` must be provided.
    """
    run = get_run()
    if not get_run():
        raise RuntimeError("You must call swanlab.data.init() before using finish()")
    if not run.is_running:
        return swanlog.error("After experiment is finished, you can't call finish() again.")
    run.finish(state, error)


def _login_in_init() -> LoginInfo:
    """在init函数中登录"""
    # 1. 如果没有登录，提示登录
    # 2. 如果登录了，发起请求，如果请求失败，重新登录，返回步骤1
    key = None
    try:
        key = get_key(os.path.join(get_swanlab_folder(), ".netrc"), get_host_api())[2]
    except KeyFileError:
        fd = sys.stdin.fileno()
        # 不是标准终端，且非jupyter环境，无法控制其回显
        if not os.isatty(fd) and not in_jupyter():
            raise KeyFileError("The key file is not found, call `swanlab.login()` or use `swanlab login` ")
    return terminal_login(key)


def _init_logdir(logdir: str) -> str:
    """
    处理通过init传入的logdir存在的一些情况
    """
    # 如果传入了logdir，则将logdir设置为环境变量，代表日志文件存放的路径
    if logdir is not None:
        try:
            logdir = check_dir_and_create(logdir)
        except ValueError:
            raise ValueError("logdir must be a str.")
        except IOError:
            raise IOError("logdir must be a path and have Write permission.")
        os.environ[ROOT] = logdir
    # 如果没有传入logdir，则使用默认的logdir, 即当前工作目录下的swanlog文件夹，但是需要保证目录存在
    else:
        logdir = os.path.abspath("swanlog")
        try:
            os.makedirs(logdir, exist_ok=True)
            if not os.access(logdir, os.W_OK):
                raise IOError
        except IOError:
            raise IOError("logdir must have Write permission.")
    return logdir


def _init_config(config: Union[dict, str]):
    """初始化传入的config参数"""
    if isinstance(config, dict) or config is None:
        return config
    swanlog.info("The parameter config is loaded from the configuration file: {}".format(config))
    return check_load_json_yaml(config, "config")


def _load_data(load_data: dict, key: str, value):
    """从load_data中加载数据，如果value不是None"""
    if value is not None:
        # tip = "The parameter {} is loaded from the configuration file: {}".format(FONT.bold(key), value)
        # print(FONT.swanlab(tip))
        return value
    d = load_data.get(key, None)
    # if d is not None:
    #     tip = "The parameter {} is loaded from the configuration file: {}".format(FONT.bold(key), d)
    #     print(FONT.swanlab(tip))
    return d
