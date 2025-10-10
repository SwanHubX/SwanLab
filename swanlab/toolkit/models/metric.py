"""
@author: cunyue
@file: metric.py
@time: 2025/7/20 16:25
@description: 指标模型
"""

import os
from typing import Optional, Dict, List, Literal, Tuple
from urllib.parse import quote

from .data import ChartType, ParseErrorInfo, MediaBuffer, ChartReference

__all__ = [
    "ColumnInfo",
    "MetricInfo",
    "MetricErrorInfo",
    "ColumnClass",
    "SectionType",
    "ColumnConfig",
    "YRange",
]


ColumnClass = Literal["CUSTOM", "SYSTEM"]
SectionType = Literal["PINNED", "HIDDEN", "PUBLIC", "SYSTEM", "CUSTOM"]
YRange = Optional[Tuple[Optional[float], Optional[float]]]


class ColumnConfig:
    """
    列信息配置
    """

    def __init__(
        self,
        y_range: YRange = None,
        chart_name: Optional[str] = None,
        chart_index: Optional[str] = None,
        metric_name: Optional[str] = None,
        metric_color: Optional[Tuple[str, str]] = None,
    ):
        """
        生成的列信息配置对象
        :param y_range: y轴范围
        :param chart_name: 图表名称
        :param chart_index: 图表索引
        :param metric_name: 指标名称
        :param metric_color: 指标颜色
        """
        self.y_range: YRange = y_range
        self.chart_name: Optional[str] = chart_name
        self.chart_index: Optional[str] = chart_index
        self.metric_name: Optional[str] = metric_name
        self.metric_color: Optional[Tuple[str, str]] = metric_color

    def clone(
        self,
        y_range: YRange = None,
        chart_name: Optional[str] = None,
        chart_index: Optional[str] = None,
        metric_name: Optional[str] = None,
        metric_color: Optional[Tuple[str, str]] = None,
    ):
        """
        克隆一个新的ColumnConfig对象，并且可以修改其中的参数
        :param y_range: y轴范围
        :param chart_name: 图表名称
        :param chart_index: 图表索引
        :param metric_name: 指标名称
        :param metric_color: 指标颜色
        :return: 新的ColumnConfig对象
        """
        return ColumnConfig(
            y_range=y_range if y_range is not None else self.y_range,
            chart_name=chart_name if chart_name is not None else self.chart_name,
            metric_name=metric_name if metric_name is not None else self.metric_name,
            chart_index=chart_index if chart_index is not None else self.chart_index,
            metric_color=metric_color if metric_color is not None else self.metric_color,
        )


class ColumnInfo:
    """
    列信息，当创建列时，会生成这个对象
    """

    def __init__(
        self,
        key: str,
        kid: str,
        name: str,
        cls: ColumnClass,
        chart_type: ChartType,
        chart_reference: ChartReference,
        section_name: Optional[str],
        section_type: SectionType,
        section_sort: Optional[int] = None,
        error: Optional[ParseErrorInfo] = None,
        config: Optional[ColumnConfig] = None,
    ):
        """
        生成的列信息对象
        :param key: 生成的列名称，作为索引键值
        :param kid: 当前实验下，列的唯一id，与保存路径等信息有关，与云端请求无关
        :param name: 列的别名
        :param cls: 列的类型，CUSTOM为自定义列，SYSTEM为系统生成列
        :param chart_type: 列对应的图表类型
        :param chart_reference: 这个列对应图表的参考系，step为步数，time为时间
        :param section_name: 列的组名
        :param section_sort: 列在section中的参考排序，不代表实际排序（目前为冗余设计）
        :param error: 列的类型错误信息
        :param config: 列的额外配置信息
        """
        self.key = key
        self.kid = kid
        self.name = name
        self.cls = cls

        self.chart_type = chart_type
        self.chart_reference = chart_reference

        self.section_name = section_name
        self.section_sort = section_sort
        self.section_type = section_type

        self.error = error
        self.config = config

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

    @property
    def key_encode(self):
        """
        对key进行url编码
        """
        return quote(self.key, safe="")


class MetricInfo:
    """
    指标信息，当新的指标数据被log时，会生成这个对象
    """

    __SUMMARY_NAME = "_summary.json"

    def __init__(
        self,
        column_info: ColumnInfo,
        metric: Optional[Dict],
        metric_buffers: Optional[List[MediaBuffer]],
        metric_summary: Optional[Dict],
        metric_step: Optional[int],
        metric_epoch: Optional[int],
        metric_file_name: Optional[str],
        swanlab_logdir: Optional[str],
        swanlab_media_dir: Optional[str],
        error: Optional[ParseErrorInfo] = None,
    ):
        """
        生成的指标信息对象
        :param column_info: 此指标对应的列信息
        :param metric: 此指标的数据
        :param metric_buffers: 此指标的媒体数据，如果为None，表示没有媒体数据
        :param metric_summary: 此指标的摘要信息
        :param metric_step: 此指标的步数
        :param metric_epoch: 此指标对应本地的行数
        :param metric_file_name: 此指标的文件名
        :param swanlab_logdir: swanlab在本次实验的log文件夹路径
        :param swanlab_media_dir: swanlab在本次实验的media文件夹路径
        :param error: 创建此指标时的错误信息
        """
        self.error = error
        self.column_info = column_info
        self.metric = metric
        self.metric_buffers = metric_buffers
        self.metric_summary = metric_summary
        self.metric_step = metric_step
        self.metric_epoch = metric_epoch
        _id = self.column_info.kid
        self.metric_file_path = None if self.is_error else os.path.join(swanlab_logdir, _id, metric_file_name)
        self.summary_file_path = None if self.is_error else os.path.join(swanlab_logdir, _id, self.__SUMMARY_NAME)
        self.swanlab_media_dir = swanlab_media_dir
        # 写入文件名称，对应上传时的文件名称：{key}/{文件名称}，文件夹名称为key
        if self.metric_buffers is not None:
            for i, buffer in enumerate(self.metric_buffers):
                buffer.file_name = "{}/{}".format(self.column_info.key_encode, metric["data"][i])

    @property
    def is_error(self) -> bool:
        """
        这条指标信息是否有错误，错误分几种：
            1. 列错误，列一开始就出现问题
            2. 重复错误
            3. 指标错误
        """
        return self.error is not None or self.column_error is not None

    @property
    def column_error(self) -> Optional[ParseErrorInfo]:
        """
        列错误信息
        """
        return self.column_info.error

    @property
    def data(self):
        """
        指标数据的data字段
        """
        if self.is_error:
            return None
        return self.metric["data"]

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self.data)


class MetricErrorInfo(MetricInfo):
    def __init__(self, column_info: ColumnInfo, error: ParseErrorInfo):
        """
        错误的指标信息，简化输入参数
        :param column_info: 此指标对应的列信息
        :param error: 创建此指标时的错误信息
        """
        super().__init__(
            column_info,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            error,
        )
