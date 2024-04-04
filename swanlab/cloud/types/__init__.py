#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/4 17:47
@File: __init__.py.py
@IDE: pycharm
@Description:
    日志类型模块，在此处定义所有的日志类型及其处理方法
"""
from enum import Enum
from watchdog.utils.dirsnapshot import DirectorySnapshotDiff
from watchdog.events import FileSystemEvent
from swanlab.log import swanlog
from typing import Callable, Tuple, Union


class LogPackageSeed:
    """
    日志包装种子，列举不同的日志包装方式、处理方式等，最终返回的应该是一个函数
    """
    ParamsType = Tuple[str, Union[FileSystemEvent]]

    @staticmethod
    def total_metadata(params: ParamsType) -> Callable:
        """
        将文件夹下的所有文件元数据上传，特指files文件夹
        :param params: 传入参数，月定为监听根路径和文件事件
        """

        async def _():
            swanlog.info(f"Upload the whole file, path: {params[0]}")

        return _

    @staticmethod
    def upload_all_metadata(params: ParamsType) -> Callable:
        """
        将文件的元数据上传
        :param params: 传入参数，月定为监听根路径和文件事件(事件可能为None）
        """

        async def _():
            swanlog.warning(f"Upload the metadata of file: {params[0]}")

        return _

    @staticmethod
    def diff_upload(params: ParamsType) -> Callable:
        """
        将文件的差异部分上传
        :param params: 传入参数，月定为监听根路径和文件事件(事件可能为None）
        """

        async def _():
            swanlog.error(f"Upload the diff of file: {FileSystemEvent}")

        return _


class LogPackageType(Enum):
    """
    日志类型枚举，涉及不同的日志包装方式
    """

    class PACKAGE:
        """日志包装类，私有类"""

        def __init__(self, on, before, after):
            self.__on = on
            self.__before = before
            self.__after = after

        @property
        def on(self):
            return self.__on

        @property
        def before(self):
            return self.__before

        @property
        def after(self):
            return self.__after

    # 媒体类型数据，只上传一次，不再管他们
    MEDIA = PACKAGE(LogPackageSeed.total_metadata, None, None)

    # 元数据类型数据，在监听之前设置钩子，监听之后全部上传
    META = PACKAGE(LogPackageSeed.total_metadata, LogPackageSeed.upload_all_metadata, None)

    # 指标类型数据，特指logs文件夹内数据，监听差异
    METRIC = PACKAGE(LogPackageSeed.total_metadata, LogPackageSeed.diff_upload, None)

    # 日志类型数据，特指console文件夹内数据，监听差异
    LOG = PACKAGE(LogPackageSeed.total_metadata, LogPackageSeed.diff_upload, None)
