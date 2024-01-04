#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:45:15
@File: swanlab/database/settings.py
@IDE: vscode
@Description:
    数据收集部分配置，此为运行时生成的配置，
"""
import os
from ..env import get_runtime_project, get_swanlog_dir, get_runtime_root


class SwanDataSettings:
    def __init__(self, exp_name: str) -> None:
        """实验名称

        Parameters
        ----------
        exp_name : str
            实验名称，实验名称应该唯一，由0-9，a-z，A-Z，" ","_","-","/"组成
            但此处不做限制
        """
        self.__exp_name: str = exp_name
        self.__root_dir: str = get_runtime_root()
        self.__proejct_path: str = get_runtime_project()

    @property
    def exp_name(self) -> str:
        """实验名称"""
        return self.__exp_name

    @property
    def root_dir(self) -> str:
        """根目录"""
        return self.__root_dir

    @property
    def project_path(self) -> str:
        """project.json文件路径"""
        return self.__proejct_path

    @property
    def output_path(self) -> str:
        """输出文件路径"""
        return os.path.join(self.exp_dir, "output.log")

    @property
    def chart_path(self) -> str:
        """图表配置路径"""
        return os.path.join(self.exp_dir, "chart.json")

    @property
    def error_path(self) -> str:
        """错误文件路径"""
        return os.path.join(self.console_dir, "error.log")

    @property
    def exp_dir(self) -> str:
        """实验文件夹路径"""
        return os.path.join(get_swanlog_dir(), self.exp_name)

    @property
    def log_dir(self) -> str:
        """记录用户日志文件夹路径"""
        return os.path.join(self.exp_dir, "logs")

    @property
    def console_dir(self) -> str:
        """记录终端日志路径"""
        return os.path.join(self.exp_dir, "console")

    @property
    def static_dir(self) -> str:
        """静态资源路径"""
        return os.path.join(self.exp_dir, "static")
