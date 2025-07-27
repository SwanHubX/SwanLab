"""
@author: cunyue
@file: key.py
@time: 2025/7/2 11:12
@description: 每一个指标都需要一个唯一的key来标识，我们将它实现为一个对象
"""

import json
import math
from typing import Optional, Tuple

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
    ChartType,
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
        self.column_info: Optional[ColumnInfo] = None
        self._media_dir = media_dir
        self._log_dir = log_dir
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
        result = data.parse()
        # 2. 解析失败则报错
        if data.error is not None:
            swanlog.warning(
                f"Log failed. Reason: Data on key '{self.key}' (step {result.step}) cannot be converted ."
                f"It should be {data.error.expected}, but it is {data.error.got}, please check the data type."
            )
            return MetricErrorInfo(column_info=self.column_info, error=data.error)
        # 3. 前后图表类型不一致则报错
        now_chart = result.chart
        expected_chart = self.column_info.chart_type
        if now_chart != expected_chart:
            swanlog.error(
                f"Data type error, key: {self.key}, "
                f"data type: {now_chart.value.column_type}, expected: {expected_chart.value.column_type}."
            )
            return MetricErrorInfo(
                column_info=self.column_info,
                error=ParseErrorInfo(
                    expected=expected_chart.value.column_type,
                    got=now_chart.value.column_type,
                    chart=result.chart,
                ),
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
        assert self.column_info is None, "Cannot create column info after creating it"
        result = data.parse()
        # SYSTEM PINNED HIDDEN 三个类型的 section 不应该传递名称
        # PUBLIC 可选是否传递名称，如果 key 包含斜杠，则使用斜杠前的部分作为section的名称
        # CUSTOM 时如果 key 包含斜杠，则使用斜杠前的部分作为section的名称，并且将 section_type 设置为 PUBLIC
        if section_type in ["PUBLIC", "CUSTOM"]:
            split_key = key.split("/")
            if len(split_key) > 1 and split_key[0]:
                # 如果key包含斜杠，则使用斜杠前的部分作为section的名称
                result.section = split_key[0]
                section_type: SectionType = "PUBLIC"
        else:
            result.section = None
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

    @classmethod
    def mock_from_remote(
        cls,
        key: str,
        column_type: str,
        column_class: ColumnClass,
        error: Optional[dict],
        media_dir: str,
        log_dir: str,
        kid: int,
        step: Optional[int],
        name: str = None,
    ) -> Tuple['SwanLabKey', ColumnInfo]:
        """
        从远程数据创建一个SwanLabKey对象，主要用于 Resume 时的图表数据标记
        mock的 column 不会被用于上传数据，只是用于标记和记录错误信息
        此对象不会用于上传数据，只是用于标记和记录错误信息
        :param key: str, key名称
        :param column_type: str, 列类型
        :param column_class: ColumnClass, 列的类型
        :param error: 列错误信息
        :param media_dir: str, 媒体目录
        :param log_dir: str, 指标目录
        :param kid: int, 代表这个列是第几个列，通常是从0开始的整数
        :param step: int, 云端最新步数，如果不传则默认为None
        :param name: str, 列的实际名称, 如果不传则默认为key
        """
        if name is None:
            name = key
        # 1. 创建对象
        key_obj = cls(key, media_dir, log_dir)
        # 2. 生成一个 ColumnInfo 对象，此 column 对象不会被用于上传
        if column_type == "FLOAT":
            chart = ChartType.LINE
        else:
            chart = getattr(ChartType, column_type, None)
        if chart is None:
            raise RuntimeError(
                f"Unknown chart type: {column_type}, maybe you need to update swanlab: pip install -U swanlab"
            )

        if error is not None:
            expected = error.get("excepted")
            got = error.get("data_class")
            if expected is None or got is None:
                raise RuntimeError(
                    f"Invalid error format: {error}, expected and got must be provided. "
                    f"Maybe you need to update swanlab: pip install -U swanlab"
                )
            error = ParseErrorInfo(expected=expected, got=got, chart=chart)
        # 3. 根据 column_class 和 chart 适配 section_type
        if column_class == "SYSTEM":
            section_type: SectionType = "SYSTEM"
        else:
            if chart is ChartType.ECHARTS:
                section_type = "CUSTOM"
            else:
                section_type = "PUBLIC"
        # 4. 创建 ColumnInfo 对象
        column_info = ColumnInfo(
            key,
            str(kid),
            name,
            column_class,
            chart,
            chart_reference="STEP",
            error=error,
            section_name=None,
            section_type=section_type,
        )
        key_obj.column_info = column_info
        # 5. 设置当前步数，resume 后不允许设置历史步数，所以需要覆盖
        if step is not None:
            for i in range(step + 1):
                key_obj.steps.add(i)
        return key_obj, column_info
