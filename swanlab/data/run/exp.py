import json
from swanlab.data.settings import SwanDataSettings
from swanlab.data.modules import BaseType, DataType
from swanlab.log import swanlog
from typing import Dict, Union, Optional
from swanlab.utils import create_time
from swanlab.utils.file import check_tag_format
from .callback import MetricInfo, ColumnInfo
from .operator import SwanLabRunOperator
from urllib.parse import quote
import os
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
        self.tags: Dict[str, SwanLabTag] = {}
        self.__operator = operator
        """回调函数"""

    def add(self, key: str, data: DataType, step: int = None) -> MetricInfo:
        """记录一条新的tag数据

        Parameters
        ----------
        key : str
            tag名称，需要检查数据类型
        data : Union[int, float, BaseType]
            tag数据，可以是浮点型，也可以是SwanLab定义的数据类型，具体与添加的图表类型有关系
            如果data是int或者float，添加图表时自动添加默认折线图，如果是BaseType，添加图表时选择对应的图表
            需要检查数据类型

        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        tag = check_tag_format(key, auto_cut=True)
        if key != tag:
            # 超过255字符，截断
            swanlog.warning(f"Tag {key} is too long, cut to 255 characters.")

        # 如果是swanlab自定义的数据类型，注入settings，方便在类内部使用
        if isinstance(data, BaseType):
            data.settings = self.settings
        # 判断tag是否存在，如果不存在则创建tag
        tag_obj: SwanLabTag = self.tags.get(tag, None)
        """
        数据库创建字段将在chart创建完成后进行
        无论格式是否正确，都会创建chart，但是如果格式不正确，不会写入日志
        """
        if tag_obj is None:
            # 将此tag对象添加到实验列表中
            tag_obj = SwanLabTag(tag, self.settings.log_dir)
            self.tags[tag] = tag_obj
            """
            由于添加图表的同时会尝试转换data的类型，但这在注入step之前
            所以此处需要手动注入依赖
            关于step，如果step格式不正确，会在tag_obj.add的时候被拦截，此处如果格式不正确直接设置为0即可
            """
            if isinstance(data, BaseType):
                data.step = 0 if step is None or not isinstance(step, int) else step
                data.tag = tag_obj.tag
            column_info = tag_obj.create_chart(tag, data)
            self.warn_type_error(tag)
            # 创建新列，生成回调
            self.__operator.on_column_create(column_info)

        # 检查tag创建时图表是否创建成功，如果失败则也没有写入数据的必要了，直接退出
        if not tag_obj.is_chart_valid:
            self.warn_chart_error(tag)
            return MetricInfo(key, tag_obj.column_info)
        key_info = tag_obj.add(data, step)
        key_info.static_dir = self.settings.static_dir
        return key_info

    def warn_type_error(self, tag: str):
        """警告类型错误
        执行此方法时需保证tag已经存在
        """
        tag_obj = self.tags[tag]
        class_name = tag_obj.now_data_type
        excepted = tag_obj.expect_data_types
        if tag_obj.is_chart_valid:
            return
        if class_name == "list":
            swanlog.error(f"Data type error, tag: {tag}, there is element of invalid data type in the list.")
        else:
            swanlog.error(f"Data type error, tag: {tag}, data type: {class_name}, excepted: {excepted}")

    def warn_chart_error(self, tag: str):
        """
        警告图表创建错误
        执行此方法时需保证tag已经存在
        """
        tag_obj = self.tags[tag]
        if tag_obj.is_chart_valid:
            return
        class_name = tag_obj.now_data_type
        if class_name == "list":
            swanlog.warning(
                f"Chart '{tag}' creation failed. "
                f"Reason: The data type in list of the tag '{tag}' is not as expected, please check the data type."
            )
        else:
            swanlog.warning(
                f"Chart '{tag}' creation failed. "
                f"Reason: The expected value type for the chart '{tag}' is one of int,"
                f"float or BaseType, but the input type is {class_name}."
            )


class SwanLabTag:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

    def __init__(self, tag: str, log_dir: str) -> None:
        """
        初始化tag对象

        Parameters
        ----------
        tag : str
            tag名称
        log_dir : str
            log文件夹路径
        """
        self.tag = tag
        self.__steps = set()
        """
        此tag已经包含的steps步骤
        """
        self.__log_dir = log_dir
        """存储文件夹路径"""
        self.data_types = [float, int]
        """默认数据类型，如果tag数据为BaseType的子类，则使用其规定的数据类型"""
        self._summary = {}
        """数据概要总结"""
        self.__data = self.__new_tags()
        """当前tag的数据"""
        self.__error = None
        """此tag在自动生成chart的时候的错误信息"""
        self.data_type = None
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

    @staticmethod
    def __is_nan(data):
        """判断data是否为nan"""
        return isinstance(data, (int, float)) and math.isnan(data)

    def add(self, data: DataType, step: int = None) -> MetricInfo:
        """添加一个数据，在内部完成数据类型转换
        如果转换失败，打印警告并退出
        并且添加数据，当前的数据保存是直接保存，后面会改成缓存形式

        Parameters
        ----------
        data : DataType
            待添加的数据
        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'

        Returns
        -------
        metric_info
            添加数据的信息，如果ok为False，则返回None
        ok
            是否添加成功，如果当前step已经存在，或者类型转换失败，返回False
        """
        if step is not None and not isinstance(step, int):
            swanlog.warning(f"Step {step} is not int, SwanLab will set it automatically.")
            step = None
        if step is None:
            step = len(self.__steps)
        if step in self.__steps:
            swanlog.warning(f"Step {step} on tag {self.tag} already exists, ignored.")
            return MetricInfo(self.tag, self.__column_info)
        more = None if not isinstance(data, BaseType) else data.get_more()
        try:
            data = self.try_convert_after_add_chart(data, step)
        except ValueError:
            swanlog.warning(
                f"Log failed. Reason: Data {data} on tag '{self.tag}' (step {step}) cannot be converted .It should be "
                f"an int, float, or a DataType, but it is {type(data)}, please check the data type."
            )
            return MetricInfo(self.tag, self.__column_info)
        is_nan = self.__is_nan(data)
        # 更新数据概要
        if not is_nan:
            if self._summary.get("max") is None or data > self._summary["max"]:
                self._summary["max"] = data
                self._summary["max_step"] = step
            if self._summary.get("min") is None or data < self._summary["min"]:
                self._summary["min"] = data
                self._summary["min_step"] = step
        self._summary["num"] = self._summary.get("num", 0) + 1
        self.__steps.add(step)
        swanlog.debug(f"Add data, tag: {self.tag}, step: {step}, data: {data}")
        if len(self.__data["data"]) >= self.__slice_size:
            self.__data = self.__new_tags()
        data = data if not is_nan else "NaN"
        new_data = self.__new_tag(step, data, more=more)
        self.__data["data"].append(new_data)
        epoch = len(self.__steps)
        mu = math.ceil(epoch / self.__slice_size)
        file_path = os.path.join(self.save_path, str(mu * self.__slice_size) + ".log")
        return MetricInfo(
            self.tag,
            self.__column_info,
            json.loads(json.dumps(new_data)),
            json.loads(json.dumps(self._summary)),
            self.data_type,
            step,
            epoch,
            metric_path=file_path,
            summary_path=os.path.join(self.save_path, "_summary.json"),
            error=False,
        )

    @property
    def save_path(self):
        """获取当前tag的保存路径

        Returns
        -------
        str
            保存路径
        """
        path = os.path.join(self.__log_dir, quote(self.tag, safe=""))
        return path

    def create_chart(self, tag: str, data: DataType) -> ColumnInfo:
        """在第一次添加tag的时候，自动创建图表和namespaces，同时写入tag和数据库信息，将创建的信息保存到数据库中
        此方法只能执行一次
        具体步骤是：
        1. 创建tag字段
        2. 如果不是baseType类型，
        :returns tag, chart_type, error

        WARNING 返回位置需与回调函数入参位置一致
        """
        if self.chart_created:
            raise ValueError(f"Chart {tag} has been created, cannot create again.")
        # 如果是非BaseType类型，写入默认命名空间，否则写入BaseType指定的命名空间
        sort = 0
        if not isinstance(data, BaseType):
            namespace, chart_type, reference, config = "default", "default", "step", None
            data_type = "default"
        else:
            namespace, types, reference, config = data.__next__()
            chart_type, self.data_types = types
            data_type = data.__class__.__name__.lower()
        # 对于namespace，如果tag的名称存在斜杠，则使用斜杠前的部分作为namespace的名称
        if "/" in tag and tag[0] != "/":
            namespace = tag.split("/")[0]
            sort = None
        """
        接下来判断tag格式的正确性，判断完毕后往source中添加一条tag记录
        在此函数中，只判断tag的格式是否正确，不记录数据
        在逻辑上只有第一次会检查tag的正确性，也就是说前端的error错误只有在第一次添加tag的时候才有可能出现
        如果第一次添加成功，后续出现错误，只会在添加的时候warning一下然后丢弃这个错误
        如果第一次添加失败，后续都不会再添加此数据，因此不会出现错误
        """
        # 如果data不是期望的data_types中的类型，尝试转换为这两个类型中的一个（优先转换为第一个）
        # 如果data是BaseType类型，会在try_convert中完成转换，此处不需要管
        error = None
        if not isinstance(data, tuple(self.data_types)):
            try:
                data = self.try_convert(data)
            except ValueError:
                # 此时代表数据异常，拿到data的__class__.__name__，生成error并保存
                if isinstance(data, BaseType):
                    class_name = data.value.__class__.__name__
                    excepted = data.expect_types()
                else:
                    class_name = data.__class__.__name__
                    excepted = [i.__name__ for i in self.data_types]
                error = {"data_class": class_name, "excepted": excepted}
        if self.__is_nan(data):
            error = {"data_class": "NaN", "excepted": [i.__name__ for i in self.data_types]}
        column_info = ColumnInfo(tag, namespace, data_type, chart_type, sort, error, reference, config)
        self.__error = error
        self.data_type = data_type
        self.__column_info = column_info
        return column_info

    @staticmethod
    def __new_tag(index, data, more: dict = None) -> dict:
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
    def __new_tags() -> dict:
        """创建一个新的tag data数据集合

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

    def try_convert(self, value: DataType):
        # 如果当前data已经是data_types中的类型，直接返回
        if isinstance(value, tuple(self.data_types)):
            return value
        if isinstance(value, BaseType):
            value = value.convert
        for data_type in self.data_types:
            try:
                if type(value) in self.data_types:
                    return value
                converted_value = data_type(value)
                return converted_value
            except (ValueError, TypeError):
                # 如果转换失败，继续尝试下一种类型
                continue
        # 如果所有类型都尝试过仍然失败，则抛出异常
        raise ValueError(f"Unable to convert {value} to any of the specified types.")

    def try_convert_after_add_chart(self, data, step):
        """尝试将data转换为data_types中的类型，如果转换失败，返回None
        如果所有类型都尝试过仍然失败，则抛出异常
        调用这个函数代表用户的图表已经创建并且传入了与之前不同的数据类型
        """
        try:
            if isinstance(data, BaseType):
                # 注入内容
                if data.step is None:
                    data.step = step
                if data.tag is None:
                    data.tag = self.tag
            return self.try_convert(data)
        except ValueError:
            raise ValueError(
                f"Unable to convert {data} to any of the specified types. This Error means that you have changed the "
                f"data type of the chart, please check the data type of the chart."
            )
