import json
from swanlab.data.settings import SwanDataSettings
from swanlab.data.modules import DataWrapper, Line
from swanlab.log import swanlog
from typing import Dict, Optional
from swanlab.utils import create_time
from .callback import MetricInfo, ColumnInfo
from .operator import SwanLabRunOperator
import math


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, settings: SwanDataSettings, operator: SwanLabRunOperator) -> None:
        """初始化实验

        Parameters
        ----------
        settings : SwanDataSettings
            全局运行时配置
        operator : SwanLabRunOperator
            操作员
        """
        self.settings = settings
        # 当前实验的所有tag数据字段
        self.keys: Dict[str, SwanLabKey] = {}
        self.__operator = operator

    def add(self, key: str, data: DataWrapper, step: int = None) -> MetricInfo:
        """记录一条新的tag数据

        Parameters
        ----------
        key : str
            key名称
        data : DataWrapper
            包装后的数据，用于数据解析

        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        # 判断tag是否存在，如果不存在则创建tag
        key_obj: SwanLabKey = self.keys.get(key, None)

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
                return MetricInfo(key, key_obj.column_info)
        data.parse(step=step, settings=self.settings, key=key)

        # ---------------------------------- 图表创建 ----------------------------------

        if key_obj is None:
            # 将此tag对象添加到实验列表中
            key_obj = SwanLabKey(key, self.settings.log_dir)
            self.keys[key] = key_obj
            # 新建图表，完成数据格式校验
            column_info = key_obj.create_chart(key, data)
            self.warn_type_error(key)
            # 创建新列，生成回调
            self.__operator.on_column_create(column_info)

        # 检查tag创建时图表是否创建成功，如果失败则也没有写入数据的必要了，直接退出
        if not key_obj.is_chart_valid:
            self.warn_chart_error(key)
            return MetricInfo(key, key_obj.column_info)
        key_info = key_obj.add(data)
        key_info.raw = data.parse().raw
        key_info.static_dir = self.settings.static_dir
        return key_info

    def warn_type_error(self, key: str):
        """警告类型错误
        执行此方法时需保证tag已经存在
        """
        tag_obj = self.keys[key]
        class_name = tag_obj.now_data_type
        excepted = tag_obj.expect_data_types
        if tag_obj.is_chart_valid:
            return
        if class_name == "list":
            swanlog.error(f"Data type error, key: {key}, there is element of invalid data type in the list.")
        else:
            swanlog.error(f"Data type error, key: {key}, data type: {class_name}, excepted: {excepted}")

    def warn_chart_error(self, key: str):
        """
        警告图表创建错误
        执行此方法时需保证key已经存在
        """
        tag_obj = self.keys[key]
        if tag_obj.is_chart_valid:
            return
        class_name = tag_obj.now_data_type
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

    def __init__(self, key: str, log_dir: str) -> None:
        """
        初始化tag对象

        Parameters
        ----------
        key : str
            列名称
        log_dir : str
            log文件夹路径
        """
        self.key = key
        self.__steps = set()
        """
        此tag已经包含的steps步骤
        """
        self.__log_dir = log_dir
        """swanlab 存储文件夹路径"""
        self.__summary = {}
        """数据概要总结"""
        self.__collection = self.__new_metric_collection()
        """当前tag的数据"""
        self.__error = None
        """此tag在自动生成chart的时候的错误信息"""
        self.chart = None
        """当前tag的数据类型，如果是BaseType类型，则为BaseType的小写类名，否则为default"""
        self.__column_info = None

    @property
    def now_data_type(self) -> Optional[str]:
        """
        当is_chart_valid为False时，返回当前数据的类型
        这指的是用户传入的数据类型，而不是转换后的数据类型
        """
        if self.__error is not None:
            return self.__error.get("data_class")
        return None

    @property
    def expect_data_types(self) -> Optional[list]:
        """
        当is_chart_valid为False时，返回期望的数据类型
        如果为True，则返回None
        """
        if self.__error is not None:
            return self.__error.get("excepted")
        return None

    @property
    def sum(self):
        """当前tag的数据总数"""
        return len(self.__steps)

    @property
    def column_info(self):
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
        """判断当前tag对应的自动创建图表是否成功
        成功则返回False，失败则返回True
        为True则一切正常
        为False则tag对应的路径不存在
        """
        return self.__error is None

    def add(self, data: DataWrapper) -> MetricInfo:
        """添加一个数据，在内部完成数据类型转换
        如果转换失败，打印警告并退出
        并且添加数据，当前的数据保存是直接保存，后面会改成缓存形式

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
        result = data.parse()
        if data.error is not None:
            swanlog.warning(
                f"Log failed. Reason: Data on key '{self.key}' (step {result.step}) cannot be converted ."
                f"It should be {data.error.expected}, but it is {data.error.got}, please check the data type."
            )
            return MetricInfo(self.key, self.__column_info)

        # 如果为Line且为NaN或者INF，不更新summary
        if not data.type == Line or result.data not in [Line.nan, Line.inf]:
            if self.__summary.get("max") is None or result.data > self.__summary["max"]:
                self.__summary["max"] = result.data
                self.__summary["max_step"] = result.step
            if self.__summary.get("min") is None or result.data < self.__summary["min"]:
                self.__summary["min"] = result.data
                self.__summary["min_step"] = result.step
        self.__summary["num"] = self.__summary.get("num", 0) + 1
        self.__steps.add(result.step)
        swanlog.debug(f"Add data, key: {self.key}, step: {result.step}, data: {result.data}")
        if len(self.__collection["data"]) >= self.__slice_size:
            self.__collection = self.__new_metric_collection()

        new_data = self.__new_metric(result.step, result.data, more=result.more)
        self.__collection["data"].append(new_data)
        epoch = len(self.__steps)
        mu = math.ceil(epoch / self.__slice_size)
        return MetricInfo(
            self.key,
            epoch=epoch,
            step=result.step,
            logdir=self.__log_dir,
            column_info=self.__column_info,
            metric=json.loads(json.dumps(new_data)),
            summary=json.loads(json.dumps(self.__summary)),
            metric_file_name=str(mu * self.__slice_size) + ".log",
        )

    def create_chart(self, key: str, data: DataWrapper) -> ColumnInfo:
        """在第一次添加tag的时候，自动创建图表和namespaces，同时写入tag和数据库信息，将创建的信息保存到数据库中
        此方法只能执行一次
        具体步骤是：
        1. 创建key字段
        2. 如果不是baseType类型，
        """
        if self.chart_created:
            raise ValueError(f"Chart {key} has been created, cannot create again.")
        result = data.parse()

        # 对于namespace，如果tag的名称存在斜杠，则使用斜杠前的部分作为namespace的名称
        if "/" in key and key[0] != "/":
            result.section = key.split("/")[0]
        # 如果出现错误
        error = None
        if data.error is not None:
            error = {"data_class": data.error.got, "excepted": data.error.expected}

        column_info = ColumnInfo(
            key=key,
            namespace=result.section,
            chart=result.chart,
            error=error,
            reference=result.reference,
            config=result.config,
        )
        self.__error = error
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
