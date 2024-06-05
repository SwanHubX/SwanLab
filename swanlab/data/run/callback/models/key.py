#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/5 16:32
@File: key.py
@IDE: pycharm
@Description:
    与Key相关的回调函数触发时的模型
"""
from typing import Union, Optional, Dict, List
from swanlab.data.modules import ChartType, WrapperErrorInfo, MediaBuffer
from urllib.parse import quote
import os


class ColumnInfo:
    """
    列信息，当创建列时，会生成这个对象
    """

    def __init__(
        self,
        key: str,
        namespace: str,
        chart: ChartType,
        sort: Optional[int] = None,
        error: Optional[WrapperErrorInfo] = None,
        reference: Optional[str] = None,
        config: Optional[Dict] = None,
    ):
        self.key = key
        """
        列的key名称
        """
        self.namespace = namespace
        """
        列的命名空间
        """
        self.chart = chart
        """
        列的图表类型
        """
        self.error = error
        """
        列的类型错误信息
        """
        self.reference = reference
        """
        列的参考对象
        """
        self.sort = sort
        """
        列在namespace中的排序
        """
        self.config = config if config is not None else {}
        """
        列的额外配置信息
        """

    @property
    def got(self):
        """
        传入的错误类型，如果列出错，返回错误类型，如果没出错，`暂时`返回None
        """
        if self.error is None:
            return None
        return self.error.got

    @property
    def expected(self):
        """
        期望的类型，如果列出错，返回期望的类型，如果没出错，`暂时`返回None
        """
        if self.error is None:
            return None
        return self.error.expected


class MetricInfo:
    """
    指标信息，当新的指标被log时，会生成这个对象
    """
    __SUMMARY_NAME = "_summary.json"

    def __init__(
        self,
        key: str,
        column_info: ColumnInfo,
        error: Optional[WrapperErrorInfo],
        metric: Union[Dict, None] = None,
        summary: Union[Dict, None] = None,
        step: int = None,
        epoch: int = None,
        logdir: str = None,
        metric_file_name: str = None,
        media_dir: str = None,
        buffers: List[MediaBuffer] = None,
    ):
        self.__error = error

        self.key = quote(key, safe="")
        """
        指标的key名称，被quote编码
        """
        self.column_info = column_info
        """
        指标对应的列信息
        """
        self.metric = metric
        """
        指标信息，error时为None
        """
        self.summary = summary
        """
        摘要信息，error时为None
        """
        self.step = step
        """
        当前指标的步数，error时为None
        """
        self.epoch = epoch
        """
        当前指标对应本地的行数，error时为None
        """
        self.metric_path = None if self.error else os.path.join(logdir, self.key, metric_file_name)
        """
        指标文件的路径，error时为None
        """
        self.summary_path = None if self.error else os.path.join(logdir, self.key, self.__SUMMARY_NAME)
        """
        摘要文件的路径，error时为None
        """
        self.media_dir = media_dir
        """
        静态文件的根文件夹
        """
        self.buffers = buffers
        """
        需要上传的媒体数据，比特流，error时为None，如果上传为非媒体类型（或Text类型），也为None
        """
        # 写入文件名称，对应上传时的文件名称
        if self.buffers is not None:
            for i, buffer in enumerate(self.buffers):
                buffer.file_name = "{}/{}".format(self.key, metric["data"][i])

    @property
    def error(self) -> bool:
        """
        这条指标信息是否有错误，错误分几种：
            1. 列错误，列一开始就出现问题
            2. 重复错误
            3. 指标错误
        """
        return self.error_info is not None or self.column_error_info is not None

    @property
    def column_error_info(self) -> Optional[WrapperErrorInfo]:
        """
        列错误信息
        """
        return self.column_info.error

    @property
    def error_info(self) -> Optional[WrapperErrorInfo]:
        """
        指标错误信息
        """
        return self.__error

    @property
    def duplicated_error(self) -> bool:
        """
        是否是重复的指标
        """
        return self.__error and self.__error.duplicated

    @property
    def data(self) -> Union[Dict, None]:
        """
        指标数据的data字段
        """
        if self.error:
            return None
        return self.metric["data"]
