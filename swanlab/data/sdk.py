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
import time
from datetime import datetime
from typing import Union, Dict, Literal, List

import platformdirs

from swanlab.env import SwanLabEnv
from swanlab.formatter import (
    check_load_json_yaml,
    check_callback_format,
    check_exp_name_format,
    check_proj_name_format,
    check_desc_format,
    check_tags_format,
    check_run_id_format,
)
from swanlab.log import swanlog
from swanlab.swanlab_settings import Settings, get_settings, set_settings
from swanlab.toolkit import SwanKitCallback
from .modules import DataType
from .run import (
    SwanLabRunState,
    SwanLabRun,
    get_run,
    get_metadata,
)
from .store import get_run_store
from .utils import (
    _init_config,
    _load_from_dict,
    _load_from_env,
    _load_list_from_env,
    _create_operator,
    _init_mode,
    should_call_after_init,
    should_call_before_init,
)
from ..core_python import create_client, auth
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
    login_info = auth.code_login(api_key, save) if api_key else auth.create_login_info(save)
    create_client(login_info)
    if api_key:
        os.environ[SwanLabEnv.API_KEY.value] = api_key


MODES = Literal["disabled", "cloud", "local", "offline"]


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
        tags: List[str] = None,
        config: Union[dict, str] = None,
        logdir: str = None,
        mode: MODES = None,
        load: str = None,
        public: bool = None,
        callbacks: List[SwanKitCallback] = None,
        settings: Settings = None,
        id: str = None,
        resume: Union[Literal['must', 'allow', 'never'], bool] = None,
        reinit: bool = None,
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
        tags : List[str], optional
            The tags of the experiment, used for labeling the current experiment.
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
            Allowed values are 'cloud', 'local', 'disabled', 'offline'.
            If the value is 'cloud', data will be uploaded to the cloud and the local log will be saved.
            If the value is 'local', data will only be saved locally and will not be uploaded to the cloud.
            If the value is 'disabled', data will not be saved or uploaded, just parsing the data.
            If the value is 'offline', data will be saved locally without uploading to the cloud.
            BTW, 'online' is an alias for 'cloud'.
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
        settings: swanlab.swanlab_settings.Settings, optional
            The settings for the current experiment.
        resume : Literal['must', 'allow', 'never'], optional
            Resume the previous run or not:
                - must: You must pass the `id` parameter and the run must exist.
                - allow: If the run exists, it will be resumed, otherwise a new run will be created.
                - never: You cannot pass the `id` parameter, and a new run will be created.
            You can also pass a boolean value, where `True` is equivalent to 'allow' and `False` is equivalent to 'never'.
            [Note that] This parameter is only valid when mode='cloud'
        id : str, optional
            The run ID of the previous run, which is used to resume the previous run.
        reinit : bool, optional
            Whether to reinitialize SwanLabRun, the default is False.
        """
        # 一个进程同时只能有一个实验在运行
        if SwanLabRun.is_started():
            run = get_run()
            if reinit:
                run.finish()
            else:
                swanlog.warning("You have already initialized a run, the init function will be ignored")
                return run
        # 注册settings
        merge_settings(settings)
        user_settings = get_settings()
        swanlog.level = kwargs.get("log_level", "info")
        # ---------------------------------- 一些变量、格式检查 ----------------------------------
        # 1. 加载参数
        if callbacks is None:
            callbacks = []
        if resume is True:
            resume = 'allow'
        elif resume is False:
            resume = 'never'
        # for https://github.com/SwanHubX/SwanLab/issues/809
        if experiment_name is None and kwargs.get("name", None) is not None:
            experiment_name = kwargs.get("name")
        if description is None and kwargs.get("notes", None) is not None:
            description = kwargs.get("notes")

        # 1.1 从文件中加载数据
        if load:
            load_data = check_load_json_yaml(load, load)
            experiment_name = _load_from_dict(load_data, "experiment_name", experiment_name)
            description = _load_from_dict(load_data, "description", description)
            tags = _load_from_dict(load_data, "tags", tags)
            config = _load_from_dict(load_data, "config", config)
            logdir = _load_from_dict(load_data, "logdir", logdir)
            mode = _load_from_dict(load_data, "mode", mode)
            project = _load_from_dict(load_data, "project", project)
            workspace = _load_from_dict(load_data, "workspace", workspace)
            public = _load_from_dict(load_data, "private", public)
            id = _load_from_dict(load_data, "id", id)
            resume = _load_from_dict(load_data, "resume", resume)
        # 1.2 初始化confi参数
        config = _init_config(config)
        # 如果config是以下几个类别之一，则抛出异常
        if isinstance(config, (int, float, str, bool, list, tuple, set)):
            raise TypeError(
                f"config: {config} (type: {type(config)}) is not a json serialized dict "
                f"(Support type is dict, MutableMapping, omegaconf.DictConfig, Argparse.Namespace and other dict-like objects), "
                "please check it"
            )
        # 1.3 从环境变量中加载参数
        workspace = _load_from_env(SwanLabEnv.WORKSPACE.value, workspace)
        project = _load_from_env(SwanLabEnv.PROJ_NAME.value, project)
        experiment_name = _load_from_env(SwanLabEnv.EXP_NAME.value, experiment_name)
        description = _load_from_env(SwanLabEnv.DESCRIPTION.value, description)
        tags = _load_list_from_env(SwanLabEnv.TAGS.value, tags)
        resume = _load_from_env(SwanLabEnv.RESUME.value, resume)
        id = _load_from_env(SwanLabEnv.RUN_ID.value, id)

        # 2. 格式校验
        # 2.1 校验项目名称
        # 默认实验名称为当前目录名
        project = project if project else os.path.basename(os.getcwd())
        p = check_proj_name_format(project)
        if len(p) != len(project):
            swanlog.warning(f"project name is too long, auto cut to {p}")
            project = p
        # 2.2 校验实验名称
        if experiment_name:
            e = check_exp_name_format(experiment_name)
            if experiment_name != e:
                swanlog.warning("The experiment name has been truncated automatically.")
                experiment_name = e
        # 2.3 校验实验描述
        if description:
            d = check_desc_format(description)
            if description != d:
                swanlog.warning("The description has been truncated automatically.")
                description = d
        # 4. 校验标签
        if tags:
            new_tags = check_tags_format(tags)
            for i in range(len(tags)):
                if tags[i] != new_tags[i]:
                    swanlog.warning("The tag has been truncated automatically.")
                    tags[i] = new_tags[i]
        # 6. 校验回调函数
        callbacks = check_callback_format(self.cbs + callbacks)
        self.cbs = []
        # 7. 校验mode参数并适配 backup 模式
        mode = "cloud" if mode == "online" else mode
        mode, login_info = _init_mode(mode)
        if mode in ["offline", "local"] and user_settings.backup is False:
            raise RuntimeError("You can't use offline or local mode with backup=False!")
        # 8. 校验 resume 与 id
        resume = resume or 'never'
        if resume == 'never':
            # 不允许传递 id
            if id is not None:
                raise RuntimeError("You can't pass id when resume=never or resume=False.")
        elif resume == 'must' or resume == 'allow':
            # 只允许在 cloud 模式下使用 resume
            if mode != "cloud":
                raise RuntimeError("You can only use resume in cloud mode.")
        if resume == "must" and id is None:
            raise ValueError('You must pass id when resume=must.')
        if id is not None:
            check_run_id_format(id)
        run_store = get_run_store()
        # ---------------------------------- 初始化swanlog文件夹 ----------------------------------
        env_key = SwanLabEnv.SWANLOG_FOLDER.value
        if user_settings.backup is False or mode == "disabled":
            # 设置 backup 为 false 时只是给用户一个没有备份的感觉，但是因为 swanlab 架构问题，必须有地方保存
            # 此时我们将日志存在系统运行时目录
            logdir = platformdirs.user_cache_dir(ensure_exists=True, appname="swanlab", appauthor="SwanHubX")
        elif logdir is None:
            # 如果传入了logdir，则将logdir设置为环境变量，代表日志文件存放的路径
            # 如果没有传入logdir，则使用默认的logdir, 即当前工作目录下的swanlog文件夹，但是需要保证目录存在
            logdir = os.environ.get(env_key) or os.path.join(os.getcwd(), "swanlog")
        logdir = os.path.abspath(logdir)
        os.environ[env_key] = logdir
        try:
            os.makedirs(logdir, exist_ok=True)
            if not os.access(logdir, os.W_OK):
                raise IOError(f"no write permission for path: {logdir}")
        except Exception as error:
            raise IOError(f"Failed to create or access logdir: {logdir}, error: {error}")
        # 如果logdir是空的，创建.gitignore文件，写入*
        if not os.listdir(logdir):
            with open(os.path.join(logdir, ".gitignore"), "w", encoding="utf-8") as f:
                f.write("*")
        run_store.swanlog_dir = logdir
        # ---------------------------------- 设置运行时配置 ----------------------------------
        # 1. 写入运行时配置
        run_store.resume = resume
        run_store.run_id = id
        run_store.project = project
        run_store.workspace = workspace
        run_store.visibility = public
        run_store.tags = tags
        run_store.description = description
        run_store.run_name = experiment_name
        run_store.swanlog_dir = logdir
        # 2. 启动操作员，注册运行实例
        operator = _create_operator(mode, login_info, callbacks)
        operator.on_init(project, workspace, public=public, logdir=logdir)
        # init结束后应该设置了一些参数
        assert run_store.run_name is not None, "Run name must be set after initialization."
        assert run_store.run_colors is not None, "Run color must be set after initialization."
        assert run_store.run_id is not None, "Run id must be set after initialization."
        assert run_store.new is not None, "Run new status must be set after initialization."
        if run_store.resume == "never":
            assert run_store.new is True, "Run new status must be True when resume never."
        if run_store.new is True:
            # 新实验时一些参数必须为 none
            assert run_store.config is None, "run_store.config should be None when new experiment."
            assert run_store.metrics is None, "run_store.metrics should be None when new experiment."
        else:
            # 恢复实验时一些参数必须存在
            assert isinstance(run_store.config, dict), "run_store.config should be dict when resuming an experiment."
            assert isinstance(run_store.metrics, dict), "run_store.metrics should be dict when resuming an experiment."
        # ---------------------------------- 初始化运行文件夹 ----------------------------------
        run_id = run_store.run_id
        run_dir = None
        while True:
            run_dir is not None and time.sleep(1)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_name = "run-{}-{}".format(timestamp, run_id)
            run_dir = os.path.join(logdir, run_name)
            try:
                os.mkdir(run_dir)
                break
            except FileExistsError:
                pass
        assert run_dir is not None, "run_dir should not be None, please check the logdir and run_id"
        run_store.run_dir = run_dir
        assert os.path.exists(run_store.run_dir), f"Run directory {run_store.run_dir} does not exist."
        os.makedirs(run_store.media_dir, exist_ok=True)
        os.makedirs(run_store.log_dir, exist_ok=True)
        os.makedirs(run_store.file_dir, exist_ok=True)
        os.makedirs(run_store.console_dir, exist_ok=True)
        # ---------------------------------- 初始化运行实例 ----------------------------------
        # 系统信息检测
        meta, monitor_funcs = None, None
        # 新实验开启系统信息检测，旧实验暂时不开启
        # 并且只在用户设置开启了元数据收集或硬件监控时才开启
        if run_store.new is True and user_settings.metadata_collect:
            meta, monitor_funcs = get_metadata(run_store.run_dir)
        run = SwanLabRun(run_config=config, operator=operator, metadata=meta, monitor_funcs=monitor_funcs)
        return run


initializer = SwanLabInitializer()

init = initializer.init

register_callbacks = initializer.register_callbacks


@should_call_after_init("You must call swanlab.init() before using log()")
def log(
    data: Dict[str, DataType],
    step: int = None,
    print_to_console: bool = False,
):
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


@should_call_before_init("You can't call merge_settings() after swanlab.init()")
def merge_settings(new_settings: Settings):
    """
    合并用户设置到全局设置
    :param new_settings: Settings对象
    :raises TypeError: 当输入不是Settings对象时抛出
    :raises RuntimeError: 当设置已被锁定时抛出
    """
    if new_settings is None:
        return

    if not isinstance(new_settings, Settings):
        raise TypeError("Expected Settings object")

    current_settings = get_settings()
    merged_data = {**current_settings.model_dump(), **new_settings.model_dump()}
    # 更新全局设置
    set_settings(Settings.model_validate(merged_data))
