"""
@author: cunyue
@file: utils.py
@time: 2025/3/11 15:03
@description: 一些为sdk.py提供的工具函数
"""

import os
from typing import Union, Optional, List, Tuple

from swankit.callback import SwanKitCallback
from swankit.env import SwanLabMode
from swankit.log import FONT

from swanlab.api import terminal_login
from swanlab.data.callbacker.cloud import CloudRunCallback
from swanlab.data.formatter import check_proj_name_format, check_load_json_yaml
from swanlab.data.run import SwanLabRun
from swanlab.data.run.helper import SwanLabRunOperator
from swanlab.env import is_interactive, SwanLabEnv
from swanlab.error import KeyFileError
from swanlab.log import swanlog
from swanlab.package import get_key, get_host_web


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


def should_call_before_init(text):
    """
    装饰器，限制必须在实验初始化前调用
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            if SwanLabRun.is_started():
                raise RuntimeError(text)
            return func(*args, **kwargs)

        return wrapper

    return decorator


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
