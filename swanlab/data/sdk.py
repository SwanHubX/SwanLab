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
from typing import Union, Dict, Literal, List

from swankit.callback import SwanKitCallback

from swanlab.api import code_login, create_http
from swanlab.env import SwanLabEnv
from swanlab.log import swanlog
from .callbacker.cloud import CloudRunCallback
from .formatter import check_load_json_yaml, check_callback_format
from .modules import DataType
from .run import (
    SwanLabRunState,
    SwanLabRun,
    register,
    get_run,
)
from .run.helper import SwanLabRunOperator
from .utils import (
    _check_proj_name,
    should_call_after_init,
    _init_config,
    _load_data,
    _create_operator,
    should_call_before_init,
)
from ..package import HostFormatter


@should_call_before_init("After calling swanlab.login(), you can't call it again.")
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


class SwanLabInitializer:
    def __init__(self):
        self.cbs: List[SwanKitCallback] = []

    @should_call_before_init("After calling swanlab.init(), you can't call it again.")
    def register_callbacks(self, callbacks: List[SwanKitCallback]) -> None:
        self.cbs += callbacks

    def init(
        self,
        project: str = None,
        workspace: str = None,
        experiment_name: str = None,
        description: str = None,
        config: Union[dict, str] = None,
        logdir: str = None,
        mode: MODES = None,
        load: str = None,
        public: bool = None,
        callbacks: List[SwanKitCallback] = None,
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
        if callbacks is None:
            callbacks = []

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
        callbacks = check_callback_format(self.cbs + callbacks)
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


initializer = SwanLabInitializer()

init = initializer.init

register_callbacks = initializer.register_callbacks


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
