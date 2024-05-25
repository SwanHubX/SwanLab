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
from ..env import get_swanlog_dir
from typing import Tuple
from swanlab.package import get_package_version


class SwanDataSettings:
    def __init__(self, run_id: str, should_save: bool) -> None:
        """实验名称

        Parameters
        ----------
        exp_name : str
            实验名称，实验名称应该唯一，由0-9，a-z，A-Z，" ","_","-","/"组成
            但此处不做限制
        run_id : str
            实验运行id，由时间戳生成，用于区分不同实验存储目录
        """
        self.__exp_name = None
        self.__exp_colors = None
        self.__description = None
        # 日志存放目录
        self.__swanlog_dir: str = get_swanlog_dir()
        # 日志存放目录的上一级目录，默认情况下这应该是项目根目录
        self.__root_dir: str = os.path.dirname(self.__swanlog_dir)
        # 实验运行id
        self.__run_id: str = run_id
        self.__version = get_package_version()
        self.__should_save = should_save

    @property
    def should_save(self):
        """
        是否应该保存实验信息
        """
        return self.__should_save

    def mkdir(self, path: str) -> None:
        """创建目录"""
        if not os.path.exists(path) and self.should_save:
            os.makedirs(path, exist_ok=True)

    @property
    def version(self) -> str:
        return self.__version

    @property
    def exp_name(self) -> str:
        """实验名称"""
        if self.__exp_name is None:
            raise ValueError("exp_name is None before set")
        return self.__exp_name

    @exp_name.setter
    def exp_name(self, exp_name: str) -> None:
        """实验名称"""
        if self.__exp_name is not None:
            raise ValueError("exp_name can only be set once")
        self.__exp_name = exp_name

    @property
    def exp_colors(self) -> Tuple[str, str]:
        """实验颜色"""
        return self.__exp_colors

    @exp_colors.setter
    def exp_colors(self, exp_colors: Tuple[str, str]) -> None:
        """实验颜色"""
        if self.__exp_colors is not None:
            raise ValueError("exp_colors can only be set once")
        self.__exp_colors = exp_colors

    @property
    def description(self) -> str:
        """实验描述"""
        return self.__description

    @description.setter
    def description(self, description: str) -> None:
        """实验描述"""
        if self.__description is not None:
            raise ValueError("description can only be set once")
        self.__description = description

    @property
    def run_id(self) -> str:
        """实验运行id"""
        return self.__run_id

    @property
    def root_dir(self) -> str:
        """根目录"""
        return self.__root_dir

    @property
    def swanlog_dir(self) -> str:
        """swanlog目录，存储所有实验的日志目录，也是runs.swanlab数据库的存储目录"""
        return self.__swanlog_dir

    @property
    def run_dir(self) -> str:
        """实验日志、信息文件夹路径"""
        return os.path.join(get_swanlog_dir(), self.run_id)

    @property
    def error_path(self) -> str:
        """错误日志文件路径"""
        return os.path.join(self.console_dir, "error.log")

    @property
    def log_dir(self) -> str:
        """记录用户日志文件夹路径"""
        return os.path.join(self.run_dir, "logs")

    @property
    def console_dir(self) -> str:
        """记录终端日志路径"""
        return os.path.join(self.run_dir, "console")

    @property
    def static_dir(self) -> str:
        """静态资源路径"""
        path = os.path.join(self.run_dir, "media")
        self.mkdir(path)
        return path

    @property
    def files_dir(self) -> str:
        """实验配置信息路径"""
        path = os.path.join(self.run_dir, "files")
        self.mkdir(path)
        return path

    @property
    def requirements_path(self) -> str:
        """实验依赖的存储文件"""
        return os.path.join(self.files_dir, "requirements.txt")

    @property
    def metadata_path(self) -> str:
        """实验环境存储文件"""
        return os.path.join(self.files_dir, "swanlab-metadata.json")

    @property
    def config_path(self) -> str:
        return os.path.join(self.files_dir, "config.yaml")
