"""
@author: cunyue
@file: utils.py
@time: 2025/3/11 15:03
@description: 一些为sdk.py提供的工具函数
"""

import os
from typing import Union, Optional, List, Any

from rich.text import Text

from swanlab.core_python import auth
from swanlab.data.run import SwanLabRun
from swanlab.data.run.helper import SwanLabRunOperator
from swanlab.env import is_interactive, SwanLabEnv
from swanlab.error import KeyFileError
from swanlab.formatter import check_load_json_yaml
from swanlab.log import swanlog
from swanlab.package import get_key, get_host_web
from swanlab.toolkit import SwanKitCallback, SwanLabMode


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
                f"Using SwanLab to track your experiments. Please refer to",
                Text('https://docs.swanlab.cn', 'yellow'),
                "for more information.",
            )
            swanlog.info("(1) Create a SwanLab account.")
            swanlog.info("(2) Use an existing SwanLab account.")
            swanlog.info("(3) Don't visualize my results.")

            web_host = get_host_web()
            # 交互选择
            swanlog.info("Enter your choice: ")
            code = input("")
            while code not in ["1", "2", "3"]:
                swanlog.warning("Invalid choice, please enter again:")
                code = input("")
            if code == "3":
                mode = "offline"
            elif code == "2":
                swanlog.info("You chose 'Use an existing swanlab account'")
                swanlog.info("Logging into", Text(web_host, 'yellow'))
                login_info = auth.terminal_login()
            elif code == "1":
                swanlog.info("You chose 'Create a swanlab account'")
                swanlog.info("Create a SwanLab account here:", Text(web_host + "/login", 'yellow'))
                login_info = auth.terminal_login()
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


def _load_from_dict(load_data: dict, key: str, value):
    """
    从load_data中加载数据，如果value不是None，则直接返回value，如果为None，则返回load_data中的key
    """
    if value is not None:
        return value
    d = load_data.get(key, None)
    return d


def _load_from_env(key: Any, value) -> Optional[str]:
    """
    从环境变量中加载数据，如果value不是None，则直接返回value，如果为None，则返回环境变量中的key
    :param key: 环境变量中的key
    :param value: 传入的value
    :return: 环境变量中的value
    """
    if value is not None:
        return value
    env_value = os.getenv(key)
    if env_value is not None:
        os.environ[key] = env_value
        return env_value


def _load_list_from_env(key: Any, value: Optional[List[str]]) -> Optional[List[str]]:
    """
    从环境变量中加载tags，如果value不是None，则直接返回value，如果为None，则返回环境变量中的key，按逗号分隔
    处理逻辑与_load_from_env类似
    例如：SWANLAB_TAGS="tag1, tag2, tag3"
    结果返回 ['tag1', 'tag2', 'tag3']
    :param key: 环境变量中的key
    :param value: 传入的value
    :return:
    """
    if value is not None:
        return value
    env_value = os.getenv(key)
    if env_value is not None:
        os.environ[key] = env_value
        return [tag.strip() for tag in env_value.split(",") if tag.strip()]


def _create_operator(
    mode: str,
    login_info: Optional[auth.LoginInfo],
    cbs: Optional[List[SwanKitCallback]],
) -> SwanLabRunOperator:
    """
    创建SwanLabRunOperator实例
    如果mode为disabled，则返回一个空的SwanLabRunOperator实例

    :param mode: 运行模式
    :param login_info: 用户登录信息，如果输入则注入到CloudRunCallback中，允许自动登录
    :param cbs: 用户传递的回调函数列表
    :return: SwanLabRunOperator, CloudRunCallback
    """
    from swanlab.data.callbacker import DisabledCallback, CloudPyCallback, OfflineCallback

    c = []
    # 1.1. 禁用模式
    if mode == SwanLabMode.DISABLED.value:
        swanlog.warning("SwanLab run disabled, the data will not be saved or uploaded.")
        return SwanLabRunOperator([DisabledCallback()])
    # 1.2. 云端模式
    elif mode == SwanLabMode.CLOUD.value:
        # 在实例化CloudRunCallback之前，注入登录信息
        CloudPyCallback.login_info = login_info
        c.append(CloudPyCallback())
    # 1.3. 本地模式
    elif mode == SwanLabMode.LOCAL.value:
        from .callbacker.local import LocalRunCallback

        # 本地模式不保存 media，由回调同步保存
        c.append(LocalRunCallback())
    # 1.4 . 备份模式
    elif mode == SwanLabMode.OFFLINE.value:
        c.append(OfflineCallback())
    # 1.5. 其他非法模式 报错，backup 模式不需要在此处理
    # 上层已经 merge_settings , get_settings().backup 与此处是否设置 backup 功能等价
    elif mode not in SwanLabMode.list():
        raise ValueError(f"Unknown mode: {mode}, please use one of {SwanLabMode.list()}")

    # 2. 合并用户传递的回调函数并注册到 SwanLabRunOperator 中使其可被调用
    # WARNING: 因为官方回调接管了 SwanLabRun 的生命周期，所以用户传递的回调函数必须在官方回调函数之前执行，也就是排在列表前面
    callbacks = cbs + c
    return SwanLabRunOperator(callbacks)
