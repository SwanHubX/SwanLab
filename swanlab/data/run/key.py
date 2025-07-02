"""
@author: cunyue
@file: key.py
@time: 2025/7/2 11:12
@description: 每一个指标都需要一个唯一的key来标识，我们将它实现为一个对象
"""

import json
import math
from typing import Optional

from swanlab.data.modules import DataWrapper, Line
from swanlab.log import swanlog
from swanlab.toolkit import (
    MetricInfo,
    ColumnInfo,
    MetricErrorInfo,
    ColumnClass,
    SectionType,
    ColumnConfig,
    create_time,
    ParseErrorInfo,
)


class SwanLabKey:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

    def __init__(
        self,
        key: str,
        media_dir: str,
        log_dir: str,
    ) -> None:
        self.key = key
        # 当前 key 包含的 step
        self.steps = set()
        # 当前 key 对应的图表信息
        self.chart = None
        self.column_info = None
        self._media_dir = media_dir
        self._log_dir = log_dir
        # 当前列数据的类型
        self._class = None
        # 当前数据概要总结
        self._summary = {}
        # 当前key的数据集合
        self._collection = self.__new_metric_collection()

    @property
    def sum(self):
        """当前tag的数据总数"""
        return len(self.steps)

    @property
    def chart_created(self):
        """判断当前tag是否已经创建了图表"""
        return self.column_info is not None

    @property
    def is_chart_valid(self) -> bool:
        """判断当前tag对应的自动创建图表是否成功"""
        return self.column_info.error is None

    def add(self, data: DataWrapper) -> MetricInfo:
        """添加一个数据，在内部完成数据类型转换
        如果转换失败，打印警告并退出
        并且添加数据，当前的数据保存是直接保存，后面会改成缓存形式
        进入此函数之前column_info必须已经创建

        Parameters
        ----------
        data : DataType
            待添加的数据

        Returns
        -------
        metric_info
            添加数据的信息，如果ok为False，则返回None
        """
        # 1. 不允许在列创建前调用此方法
        assert self.column_info is not None, "Column info is None, please create column info first"
        assert self._class is not None, "Class is None, please create column info first"
        result = data.parse()
        # 2. 解析失败则报错
        if data.error is not None:
            swanlog.warning(
                f"Log failed. Reason: Data on key '{self.key}' (step {result.step}) cannot be converted ."
                f"It should be {data.error.expected}, but it is {data.error.got}, please check the data type."
            )
            return MetricErrorInfo(column_info=self.column_info, error=data.error)
        # 3. 前后数据类型不统一则报错
        data_type, data_expected = data.get_class(), self._class
        if data_type != data_expected:
            swanlog.error(
                f"Data type error, key: {self.key}, "
                f"data type: {data_type.__name__}, expected: {data_expected.__name__}."
            )
            return MetricErrorInfo(
                column_info=self.column_info,
                error=ParseErrorInfo(expected=data_expected.__name__, got=data_type.__name__, chart=result.chart),
            )
        # 4. 更新 summary 并添加数据
        # 如果为Line且为NaN或者INF，不更新summary
        r = result.strings or result.float
        if not data.type == Line or r not in [Line.nan, Line.inf]:
            if self._summary.get("max") is None or r > self._summary["max"]:
                self._summary["max"] = r
                self._summary["max_step"] = result.step
            if self._summary.get("min") is None or r < self._summary["min"]:
                self._summary["min"] = r
                self._summary["min_step"] = result.step
        self._summary["num"] = self._summary.get("num", 0) + 1
        self.steps.add(result.step)
        swanlog.debug(f"Add data, key: {self.key}, step: {result.step}, data: {r}")
        if len(self._collection["data"]) >= self.__slice_size:
            self._collection = self.__new_metric_collection()

        new_data = self.__new_metric(result.step, r, more=result.more)
        self._collection["data"].append(new_data)
        epoch = len(self.steps)
        mu = math.ceil(epoch / self.__slice_size)
        return MetricInfo(
            column_info=self.column_info,
            metric=json.loads(json.dumps(new_data)),
            metric_summary=json.loads(json.dumps(self._summary)),
            metric_epoch=epoch,
            metric_step=result.step,
            metric_buffers=result.buffers,
            metric_file_name=str(mu * self.__slice_size) + ".log",
            swanlab_logdir=self._log_dir,
            swanlab_media_dir=self._media_dir if result.buffers else None,
        )

    def create_column(
        self,
        key: str,
        name: Optional[str],
        column_class: ColumnClass,
        column_config: Optional[ColumnConfig],
        section_type: SectionType,
        data: DataWrapper,
        num: int,
    ) -> ColumnInfo:
        """
        创建列信息，对当前key的基本信息做一个记录
        :param key: str, key名称
        :param name: str, key的实际名称
        :param column_class: str, key的类型
        :param column_config: ColumnConfig, key的配置
        :param section_type: str, key的组类型
        :param data: DataType, 数据
        :param num: 创建此列之前的列数量
        """
        if self.chart_created:
            raise ValueError(f"Chart {key} has been created, cannot create again.")
        result = data.parse()

        # 如果section_type不为public,则section_name为None
        if section_type != "PUBLIC":
            result.section = None
        # 如果斜杠，则使用斜杠前的部分作为section的名称
        elif "/" in key and key[0] != "/":
            result.section = key.split("/")[0]

        column_info = ColumnInfo(
            key=key,
            kid=str(num),
            name=name,
            cls=column_class,
            config=column_config,
            chart_type=result.chart,
            section_name=result.section,
            section_type=section_type,
            section_sort=None,
            chart_reference=result.reference,
            error=data.error,
        )
        self.chart = result.chart.value
        self._class = data.get_class()
        self.column_info = column_info
        return column_info

    @staticmethod
    def __new_metric(index, data, more: dict = None) -> dict:
        """创建一个新的data数据，实际上是一个字典，包含一些默认信息

        Parameters
        ----------
        index : int
            步数
        data : DataType
            数据
        more : dict, optional
            更多的数据，如果有的话
        """
        if more is None:
            return {
                "index": int(index),
                "data": data,
                "create_time": create_time(),
            }
        else:
            return {
                "index": int(index),
                "data": data,
                "create_time": create_time(),
                "more": more,
            }

    @staticmethod
    def __new_metric_collection() -> dict:
        """创建一个新的key data数据集合

        Returns
        -------
        dict
            返回一个新的data数据集合
        """
        time = create_time()
        return {
            "create_time": time,
            "update_time": time,
            "data": [],
        }
