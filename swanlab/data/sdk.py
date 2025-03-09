#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 18:00:04
@File: swanlab/data/sdk.py
@IDE: vscode
@Description:
    在此处封装swanlab在日志记录模式下的各种接口
"""
import os
from typing import Optional, Union, Dict, Tuple, Literal, List

from swankit.callback import SwanKitCallback
from swankit.env import SwanLabMode
from swankit.log import FONT

from swanlab.api import code_login, terminal_login, create_http
from swanlab.env import SwanLabEnv, is_interactive
from swanlab.log import swanlog
from .callbacker.cloud import CloudRunCallback
from .formatter import check_load_json_yaml, check_proj_name_format, check_callback_format
from .modules import DataType
from .run import (
    SwanLabRunState,
    SwanLabRun,
    register,
    get_run,
)
from .run.helper import SwanLabRunOperator
from ..error import KeyFileError
from ..package import get_key, get_host_web, HostFormatter


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


def login(api_key: str = None, host: str = None, web_host: str = None, save: bool = False):
    """
    Login to SwanLab Cloud. If you already have logged in, you can use this function to relogin.
    Every time you call this function, the previous login information will be overwritten.

    [Note that] this function should be called before `init`.

    :param api_key: str, optional
        authentication key, if not provided, the key will be read from the key file.
    :param host: str, optional
        api host, if not provided, the default host will be used.
    :param web_host: str, optional
        web host, if not provided, the default host will be used.
    :param save: bool, optional
        whether to save the api key to the key file, default is False
    :return: LoginInfo
    """
    if SwanLabRun.is_started():
        raise RuntimeError("You must call swanlab.login() before using init()")
    # 检查host是否合法，并格式化，注入到环境变量中
    HostFormatter(host, web_host)()
    # 登录，初始化http对象
    login_info = code_login(api_key, save) if api_key else CloudRunCallback.create_login_info(save)
    create_http(login_info)
    if api_key:
        os.environ[SwanLabEnv.API_KEY.value] = api_key


MODES = Literal["disabled", "cloud", "local"]


def init(
    project: str = None,
    workspace: str = None,
    experiment_name: str = None,
    description: str = None,
    config: Union[dict, str] = None,
    logdir: str = None,
    mode: MODES = None,
    load: str = None,
    public: bool = None,
    callbacks: Optional[Union[SwanKitCallback, List[SwanKitCallback]]] = None,
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
    mode : str, optional
        Allowed values are 'cloud', 'cloud-only', 'local', 'disabled'.
        If the value is 'cloud', the data will be uploaded to the cloud and the local log will be saved.
        If the value is 'cloud-only', the data will only be uploaded to the cloud and the local log will not be saved.
        If the value is 'local', the data will only be saved locally and will not be uploaded to the cloud.
        If the value is 'disabled', the data will not be saved or uploaded, just parsing the data.
    load : str, optional
        If you pass this parameter,SwanLab will search for the configuration file you specified
        (which must be in JSON or YAML format)
        and automatically fill in some explicit parameters of this function for you
        (excluding parameters in `**kwargs` and the parameters if they are None).
        In terms of priority, if the parameters passed to init are `None`,
        SwanLab will attempt to replace them from the configuration file you provided;
        otherwise, it will use the parameters you passed as the definitive ones.
    public : bool, optional
        Whether the project can be seen by anyone, the default is None, which means the project is private.
        Only available in cloud mode while the first time you create the project.
    callbacks : Union[SwanKitCallback, List[SwanKitCallback]], optional
        The callback function that will be triggered when the experiment is finished.
        If you provide a list, all the callback functions in the list will be triggered in order.
    """
    if SwanLabRun.is_started():
        swanlog.warning("You have already initialized a run, the init function will be ignored")
        return get_run()
    # ---------------------------------- 一些变量、格式检查 ----------------------------------

    # for https://github.com/SwanHubX/SwanLab/issues/809
    if experiment_name is None and kwargs.get("name", None) is not None:
        experiment_name = kwargs.get("name")
    if description is None and kwargs.get("notes", None) is not None:
        description = kwargs.get("notes")

    # 从文件中加载数据
    if load:
        load_data = check_load_json_yaml(load, load)
        experiment_name = _load_data(load_data, "experiment_name", experiment_name)
        description = _load_data(load_data, "description", description)
        config = _load_data(load_data, "config", config)
        logdir = _load_data(load_data, "logdir", logdir)
        mode = _load_data(load_data, "mode", mode)
        project = _load_data(load_data, "project", project)
        workspace = _load_data(load_data, "workspace", workspace)
        public = _load_data(load_data, "private", public)
    # FIXME 没必要多一个函数
    project = _check_proj_name(project if project else os.path.basename(os.getcwd()))  # 默认实验名称为当前目录名
    callbacks = check_callback_format(callbacks)
    # ---------------------------------- 启动操作员 ----------------------------------
    operator, c = _create_operator(mode, public, callbacks)
    exp_num = SwanLabRunOperator.parse_return(
        operator.on_init(project, workspace, logdir=logdir),
        key=c.__str__() if c else None,
    )
    # 初始化confi参数
    config = _init_config(config)
    # ---------------------------------- 实例化实验 ----------------------------------
    # 注册实验
    run = register(
        project_name=project,
        experiment_name=experiment_name,
        description=description,
        run_config=config,
        log_level=kwargs.get("log_level", "info"),
        exp_num=exp_num,
        operator=operator,
    )
    return run


def should_call_after_init(text):
    """
    装饰器，限制必须在实验初始化后调用
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            if not SwanLabRun.is_started():
                raise RuntimeError(text)
            return func(*args, **kwargs)

        return wrapper

    return decorator


@should_call_after_init("You must call swanlab.init() before using log()")
def log(data: Dict[str, DataType], step: int = None, print_to_console: bool = False):
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
    print_to_console : bool, optional
        Whether to print the data to the console, the default is False.
    """
    run = get_run()
    ll = run.log(data, step)
    print_to_console and print(ll)
    return ll


@should_call_after_init("You must call swanlab.init() before using finish()")
def finish(state: SwanLabRunState = SwanLabRunState.SUCCESS, error=None):
    """
    Finish the current run and close the current experiment
    Normally, swanlab will run this function automatically,
    but you can also execute it manually and mark the experiment as 'completed'.
    Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
    If you mark the experiment as 'CRASHED' manually, `error` must be provided.
    """
    run = get_run()
    if not run.running:
        return swanlog.error("After experiment is finished, you can't call finish() again.")
    run.finish(state, error)


def _init_mode(mode: str = None):
    """
    初始化mode参数
    从环境变量中提取默认的mode参数，如果传入的mode参数不为None，则使用环境变量中的mode参数，否则使用传入的mode参数
    传入的mode必须为SwanLabMode枚举中的一个值，否则报错ValueError
    如果环境变量和传入的mode参数都为None，则默认为cloud

    从环境变量中提取mode参数以后，还有一步让用户选择运行模式的交互，详见issue： https://github.com/SwanHubX/SwanLab/issues/632

    :param mode: str, optional
        传入的mode参数
    :return: str mode
    :raise ValueError: mode参数不合法
    """
    allowed = [m.value for m in SwanLabMode]
    mode_key = SwanLabEnv.MODE.value
    mode_value = os.environ.get(mode_key)
    if mode_value is not None and mode is not None:
        swanlog.warning(f"The environment variable {mode_key} will be overwritten by the parameter mode")
    mode = mode_value if mode is None else mode
    if mode is not None and mode not in allowed:
        raise ValueError(f"`mode` must be one of {allowed}, but got {mode}")
    mode = "cloud" if mode is None else mode
    # 如果mode为cloud，且没找到 api key或者未登录，则提示用户输入
    try:
        get_key()
        no_api_key = False
    except KeyFileError:
        no_api_key = True
    login_info = None
    # 三选一只允许登录官方的host，除非在此之前手动设置了环境变量
    # 详见 https://github.com/SwanHubX/SwanLab/issues/792#issuecomment-2603959483
    if mode == "cloud" and no_api_key:
        # 判断当前进程是否在交互模式下
        if is_interactive():
            swanlog.info(
                f"Using SwanLab to track your experiments. Please refer to {FONT.yellow('https://docs.swanlab.cn')} for more information."
            )
            swanlog.info("(1) Create a SwanLab account.")
            swanlog.info("(2) Use an existing SwanLab account.")
            swanlog.info("(3) Don't visualize my results.")

            web_host = get_host_web()
            # 交互选择
            tip = FONT.swanlab("Enter your choice: ")
            code = input(tip)
            while code not in ["1", "2", "3"]:
                swanlog.warning("Invalid choice, please enter again.")
                code = input(tip)
            if code == "3":
                mode = "local"
            elif code == "2":
                swanlog.info("You chose 'Use an existing swanlab account'")
                swanlog.info("Logging into " + FONT.yellow(web_host))
                login_info = terminal_login()
            elif code == "1":
                swanlog.info("You chose 'Create a swanlab account'")
                swanlog.info("Create a SwanLab account here: " + FONT.yellow(web_host + "/login"))
                login_info = terminal_login()
            else:
                raise ValueError("Invalid choice")
        # 如果不在就不管

    os.environ[mode_key] = mode
    return mode, login_info


def _init_config(config: Union[dict, str]):
    """初始化传入的config参数"""
    if isinstance(config, str):
        swanlog.info("The parameter config is loaded from the configuration file: {}".format(config))
        return check_load_json_yaml(config, "config")

    return config


def _load_data(load_data: dict, key: str, value):
    """
    从load_data中加载数据，如果value不是None，则直接返回value，如果为None，则返回load_data中的key
    """
    if value is not None:
        return value
    d = load_data.get(key, None)
    return d


def _create_operator(
    mode: str,
    public: bool,
    cbs: Optional[Union[SwanKitCallback, List[SwanKitCallback]]],
) -> Tuple[SwanLabRunOperator, Optional[CloudRunCallback]]:
    """
    创建SwanLabRunOperator实例
    如果mode为disabled，则返回一个空的SwanLabRunOperator实例和None

    :param mode: 运行模式
    :param public: 是否公开
    :return: SwanLabRunOperator, CloudRunCallback
    """
    mode, login_info = _init_mode(mode)
    CloudRunCallback.login_info = login_info

    if mode == SwanLabMode.DISABLED.value:
        swanlog.warning("SwanLab run disabled, the data will not be saved or uploaded.")
        return SwanLabRunOperator(), None
    # 云端模式
    if mode == SwanLabMode.CLOUD.value:
        c = CloudRunCallback(public)
    # 本地模式
    else:
        from .callbacker.local import LocalRunCallback

        c = LocalRunCallback()

    callbacks = [c] + cbs
    return SwanLabRunOperator(callbacks), c
