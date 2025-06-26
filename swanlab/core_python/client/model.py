"""
@author: cunyue
@file: model.py
@time: 2025/6/16 14:55
@description: 实验、项目元信息
"""

from typing import Optional


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
