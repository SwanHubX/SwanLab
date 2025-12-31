"""
@author: cunyue
@file: utils.py
@time: 2025/12/31 13:29
@description: 客户端工具函数
"""

from typing import Tuple, Optional, Union

import requests
from urllib3.exceptions import (
    MaxRetryError,
    TimeoutError,
    NewConnectionError,
    ConnectionError,
    ReadTimeoutError,
    ConnectTimeoutError,
)

from swanlab.error import NetworkError


def safe_request(func):
    """
    在一些接口中我们不希望线程奔溃，而是返回一个错误对象
    """

    def wrapper(*args, **kwargs) -> Tuple[Optional[Union[dict, str]], Optional[Exception]]:
        try:
            # 在装饰器中调用被装饰的异步函数
            result = func(*args, **kwargs)
            return result, None
        except requests.exceptions.Timeout:
            return None, NetworkError()
        except requests.exceptions.ConnectionError:
            return None, NetworkError()
        # Catch urllib3 specific errors
        except (
            MaxRetryError,
            TimeoutError,
            NewConnectionError,
            ConnectionError,
            ReadTimeoutError,
            ConnectTimeoutError,
        ):
            return None, NetworkError()
        except Exception as e:
            return None, e

    return wrapper


class ProjectInfo:
    def __init__(self, data: dict):
        self.__data = data

    @property
    def cuid(self):
        return self.__data["cuid"]

    @property
    def name(self):
        return self.__data["name"]

    @property
    def history_exp_count(self):
        return self.__data.get('_count', {'experiments': 0})["experiments"]


class ExperimentInfo:

    def __init__(self, data: dict):
        self.__data = data

    @property
    def flag_id(self):
        """
        此实验的标志ID，标志上传时的实验会话
        """
        return self.__data["flagId"]

    @property
    def cuid(self):
        return self.__data["cuid"]

    @property
    def name(self):
        return self.__data["name"]

    @property
    def config(self) -> dict:
        """
        此实验的配置，用于 resume 时同步 config
        """
        return self.__data.get("profile", {}).get("config", {})

    @property
    def root_proj_cuid(self) -> Optional[str]:
        """
        根项目的cuid（上传实验时需要）
        """
        return self.__data.get("rootProId", None)

    @property
    def root_exp_cuid(self) -> Optional[str]:
        """
        根实验的cuid（上传实验时需要）
        """
        return self.__data.get("rootExpId", None)
