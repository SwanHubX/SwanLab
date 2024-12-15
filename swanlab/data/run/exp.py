import json
import math
from typing import Dict, Optional

from swankit.callback.models import MetricInfo, ColumnInfo, MetricErrorInfo, ColumnClass, SectionType, ColumnConfig
from swankit.core import SwanLabSharedSettings
from swankit.env import create_time

from swanlab.data.modules import DataWrapper, Line
from swanlab.log import swanlog
from .helper import SwanLabRunOperator


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, settings: SwanLabSharedSettings, operator: SwanLabRunOperator) -> None:
        """初始化实验

        Parameters
        ----------
        settings : SwanLabSharedSettings
            全局运行时配置
        operator : SwanLabRunOperator
            操作员
        """
        self.settings = settings
        # 当前实验的所有tag数据字段
        self.keys: Dict[str, SwanLabKey] = {}
        # TODO 操作员不传递给实验
        self.__operator = operator

    def __add(
        self,
        key: str,
        name: Optional[str],
        column_class: ColumnClass,
        column_config: Optional[ColumnConfig],
        section_type: SectionType,
        data: DataWrapper,
        step: int = None,
    ) -> MetricInfo:
        """记录一条新的key数据

        Parameters
        ----------
        key : str
            key的云端唯一标识
        name : str
            key的实际名称
        column_class : str
            列类型，CUSTOM为自定义key，SYSTEM为系统key
        section_type : str
            key的组类型
        data : DataWrapper
            包装后的数据，用于数据解析
        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        key_index = f"{column_class}-{key}"
        # 判断tag是否存在，如果不存在则创建tag
        key_obj: SwanLabKey = self.keys.get(key_index, None)

        # ---------------------------------- 包装器解析 ----------------------------------

        if step is not None and not isinstance(step, int):
            swanlog.warning(f"Step {step} is not int, SwanLab will set it automatically.")
            step = None
        if key_obj is None:
            step = 0 if step is None or not isinstance(step, int) else step
        else:
            step = len(key_obj.steps) if step is None else step
            if step in key_obj.steps:
                swanlog.warning(f"Step {step} on key {key} already exists, ignored.")
                return MetricErrorInfo(column_info=key_obj.column_info, error=DataWrapper.create_duplicate_error())
        data.parse(step=step, key=key)

        # ---------------------------------- 图表创建 ----------------------------------

        if key_obj is None:
            num = len(self.keys)
            # 将此tag对象添加到实验列表中
            key_obj = SwanLabKey(key, self.settings)
            self.keys[key_index] = key_obj
            # 新建图表，完成数据格式校验
            column_info = key_obj.create_column(
                key,
                name,
                column_class,
                column_config,
                section_type,
                data,
                num,
            )
            self.warn_type_error(key_index, key)
            # 创建新列，生成回调
            self.__operator.on_column_create(column_info)

        # 检查tag创建时图表是否创建成功，如果失败则也没有写入数据的必要了，直接退出
        if not key_obj.is_chart_valid:
            self.warn_chart_error(key_index, key)
            return MetricErrorInfo(key_obj.column_info, error=key_obj.column_info.error)
        key_info = key_obj.add(data)
        key_info.buffers = data.parse().buffers
        key_info.media_dir = self.settings.media_dir
        return key_info

    def add(
        self,
        data: DataWrapper,
        key: str,
        name: str = None,
        column_class: ColumnClass = 'CUSTOM',
        column_config: Optional[ColumnConfig] = None,
        section_type: SectionType = "PUBLIC",
        step: int = None,
    ) -> MetricInfo:
        """记录一条新的key数据
        Parameters
        ----------
        data : DataWrapper
            包装后的数据，用于数据解析
        key : str
            列的云端唯一标识
        name : str
            列的实际名称, 默认与key相同
        column_class : str, optional
            列的类型
        column_config : Optional[ColumnConfig], optional
            列的额外配置信息
        section_type : str, optional
            key的组类型
        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        m = self.__add(key, name, column_class, column_config, section_type, data, step)
        self.__operator.on_metric_create(m)
        return m

    def warn_type_error(self, key_index: str, key: str):
        """警告类型错误
        执行此方法时需保证key已经存在
        """
        tag_obj = self.keys[key_index]
        class_name = tag_obj.column_info.got
        expected = tag_obj.column_info.expected
        if tag_obj.is_chart_valid:
            return
        if class_name == "list":
            swanlog.error(f"Data type error, key: {key}, there is element of invalid data type in the list.")
        else:
            swanlog.error(f"Data type error, key: {key}, data type: {class_name}, expected: {expected}")

    def warn_chart_error(self, key_index: str, key: str):
        """
        警告图表创建错误
        执行此方法时需保证key已经存在
        """
        tag_obj = self.keys[key_index]
        if tag_obj.is_chart_valid:
            return
        class_name = tag_obj.column_info.got
        if class_name == "list":
            swanlog.warning(
                f"Chart '{key}' creation failed. "
                f"Reason: The data type in list of the key '{key}' is not as expected, please check the data type."
            )
        else:
            swanlog.warning(
                f"Chart '{key}' creation failed. "
                f"Reason: The expected value type for the chart '{key}' is one of int,"
                f"float or BaseType, but the input type is {class_name}."
            )


class SwanLabKey:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

    def __init__(self, key: str, settings: SwanLabSharedSettings) -> None:
        """
        初始化tag对象

        Parameters
        ----------
        key : str
            列名称
        settings : SwanLabSharedSettings
            全局运行时配置
        """
        self.key = key
        self.__steps = set()
        """
        此tag已经包含的steps步骤
        """
        self.__settings = settings
        self.__summary = {}
        """数据概要总结"""
        self.__collection = self.__new_metric_collection()
        """当前tag的数据"""
        self.chart = None
        """当前tag的数据类型，如果是BaseType类型，则为BaseType的小写类名，否则为default"""
        self.__column_info = None

    @property
    def sum(self):
        """当前tag的数据总数"""
        return len(self.__steps)

    @property
    def column_info(self) -> Optional[ColumnInfo]:
        """获取当前tag对应的ColumnInfo"""
        return self.__column_info

    @property
    def steps(self):
        """获取当前tag的所有步数"""
        return self.__steps

    @property
    def chart_created(self):
        """判断当前tag是否已经创建了图表"""
        return self.__column_info is not None

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
        ok
            是否添加成功，如果当前step已经存在，或者类型转换失败，返回False
        """
        if self.__column_info is None:
            raise ValueError("Column info is None, please create column info first")
        result = data.parse()
        if data.error is not None:
            swanlog.warning(
                f"Log failed. Reason: Data on key '{self.key}' (step {result.step}) cannot be converted ."
                f"It should be {data.error.expected}, but it is {data.error.got}, please check the data type."
            )
            return MetricErrorInfo(column_info=self.__column_info, error=data.error)

        # 如果为Line且为NaN或者INF，不更新summary
        r = result.strings or result.float

        if not data.type == Line or r not in [Line.nan, Line.inf]:
            if self.__summary.get("max") is None or r > self.__summary["max"]:
                self.__summary["max"] = r
                self.__summary["max_step"] = result.step
            if self.__summary.get("min") is None or r < self.__summary["min"]:
                self.__summary["min"] = r
                self.__summary["min_step"] = result.step
        self.__summary["num"] = self.__summary.get("num", 0) + 1
        self.__steps.add(result.step)
        swanlog.debug(f"Add data, key: {self.key}, step: {result.step}, data: {r}")
        if len(self.__collection["data"]) >= self.__slice_size:
            self.__collection = self.__new_metric_collection()

        new_data = self.__new_metric(result.step, r, more=result.more)
        self.__collection["data"].append(new_data)
        epoch = len(self.__steps)
        mu = math.ceil(epoch / self.__slice_size)
        return MetricInfo(
            column_info=self.__column_info,
            metric=json.loads(json.dumps(new_data)),
            metric_summary=json.loads(json.dumps(self.__summary)),
            metric_epoch=epoch,
            metric_step=result.step,
            metric_buffers=result.buffers,
            metric_file_name=str(mu * self.__slice_size) + ".log",
            swanlab_logdir=self.__settings.log_dir,
            swanlab_media_dir=self.__settings.media_dir if result.buffers else None,
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
        self.__column_info = column_info
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
