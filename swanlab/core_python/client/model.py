"""
@author: cunyue
@file: model.py
@time: 2025/6/16 14:55
@description: 实验、项目元信息
"""


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
    def cuid(self):
        return self.__data["cuid"]

    @property
    def name(self):
        return self.__data["name"]
