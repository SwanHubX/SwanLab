#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-28 22:12:37
@File: swanlab/database/modules/base.py
@IDE: vscode
@Description:
    SwanLab数据基类
    数据的转换流程是：
    1. 如果传入的数据在这之前已经有key，跳转步骤3
    2. 没有key的数据，需要判断当前数据是否满足要求，因此需要转换一次，跳转步骤3
    3. 转换一次，判断与之前的数据是否一致，因此还需要转换一次
    可以看到有两次转换，如果每次都是耗时操作，那么会导致性能问题，所以对于单个类而言，上次转换的结果会被保存
"""
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, ByteString, Union, Tuple
from swanlab.log import swanlog
from enum import Enum
from io import BytesIO
import hashlib
import math
import io


class DataSuite:
    """
    BaseType在处理时的通用工具套件
    """

    @staticmethod
    def get_hash_by_bytes(file_byte: ByteString, chunk: int = 4096) -> str:
        """
        计算并返回给定字节流的SHA-256哈希值
        :param file_byte: 字节流
        :param chunk: 一次读取的字节数，用于加快处理速度
        """

        hash_sha256 = hashlib.sha256()
        # 更新哈希值，为了处理速度，只读取一小块
        hash_sha256.update(file_byte[:chunk])
        return hash_sha256.hexdigest()

    @staticmethod
    def get_hash_by_ndarray(array) -> str:
        """
        计算并返回给定NumPy数组的SHA-256哈希值
        :param array: NumPy数组
        """

        hash_sha256 = hashlib.sha256()
        # 将NumPy数组转换为字节串，然后更新哈希值
        hash_sha256.update(array.tobytes())
        return hash_sha256.hexdigest()

    @staticmethod
    def get_hash_by_pil(image, f="PNG") -> str:
        """计算并返回给定PIL.Image对象的SHA-256哈希值。
        :param image: PIL.Image对象
        :param f: 图像格式，默认为PNG，可选为JPEG
        """

        hash_sha256 = hashlib.sha256()
        # 将图像转换为字节数据
        with io.BytesIO() as buffer:
            image.save(buffer, format=f)
            hash_sha256.update(buffer.getvalue())
        return hash_sha256.hexdigest()

    @staticmethod
    def is_nan(data):
        """判断data是否为nan"""
        return isinstance(data, (int, float)) and math.isnan(data)

    @staticmethod
    def is_inf(data):
        """判断data是否为inf"""
        return isinstance(data, (int, float)) and math.isinf(data)

    @staticmethod
    def check_caption(caption) -> Optional[str]:
        """
        检查caption是否符合要求
        :param caption: 标题
        :return: 返回处理后的标题
        :raises TypeError: caption不符合要求
        """
        if isinstance(caption, str):
            caption = caption
        # 如果类型是数字，则转换为字符串
        elif isinstance(caption, (int, float)):
            caption = str(caption)
        # 如果类型是None，则转换为默认字符串
        elif caption is None:
            caption = None
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption.strip() if caption else None


class DynamicProperty:
    """
    动态属性，一开始为None，动态设置
    但是可以保证在转换时都已经被设置
    """

    def __init__(self):
        # 保存的key
        self.key: Optional[str] = None
        # 保存的step
        self.step: Optional[int] = None

    def inject(self, key: str, step: int):
        """
        注入属性
        :param key: key名称
        :param step: 当前步骤
        """
        self.key = key
        self.step = step


class BaseType(ABC, DynamicProperty):
    """
    SwanLab数据基类，完成数据解析，拿到图表类型、namespace和解析后的数据
    所有的数据类都只负责解析，并返回相应的内容
    """

    class Chart(Enum):
        """
        图表选择器
        """

        class ChartItem:
            """
            图表项
            """

            def __init__(self, chart_type: str, column_type: str):
                self.chart_type = chart_type
                """
                本地标注的图表类型
                """
                self.column_type = column_type
                """
                上传列信息时，标注的图表类型
                """

        LINE = ChartItem("line", "FLOAT")

        IMAGE = ChartItem("image", "IMAGE")

        AUDIO = ChartItem("audio", "AUDIO")

        TEXT = ChartItem("text", "TEXT")

    # ---------------------------------- 需要子类实现的方法 ----------------------------------

    @abstractmethod
    def parse(self) -> Tuple[Union[str, float], Optional[BytesIO]]:
        """
        解析数据
        :raises DataTypeError: 传入的数据类型不符合要求
        """
        pass

    @abstractmethod
    def get_chart(self) -> "Chart":
        """
        获取图表类型
        """
        pass

    # ---------------------------------- 可选覆写方法 ----------------------------------
    # noinspection PyMethodMayBeStatic
    def get_section(self) -> str:
        """
        获取默认的section名称，可能在后续处理时，与实际的section名称不同
        """
        return "default"

    # noinspection PyMethodMayBeStatic
    def get_config(self) -> Optional[Dict]:
        """
        获取图表的config配置信息，应该返回一个字典，或者为None
        """
        return None

    def get_more(self) -> Optional[Dict]:
        """
        代表当前步骤的此数据支持标注的更多内容，应该返回一个字典，或者为None
        """
        return None


class MediaType(BaseType):  # noqa
    """
    媒体类型，用于区分标量和媒体，不用做实例化，应该由子类继承
    """
    pass


class MediaBuffer(BytesIO):
    def __init__(self):
        super().__init__()
        self.__file_name = None

    @property
    def file_name(self):
        if self.__file_name is None:
            raise ValueError("file_name is not set")
        return self.__file_name

    @file_name.setter
    def file_name(self, value):
        if not isinstance(value, str) or not value:
            raise TypeError(f"Expected str, but got {type(value)}")
        # if self.__file_name is not None:
        #     # 此时意味着使用类似 [] * 2 的操作复制了多个相同的实例，这允许，但不推荐
        #     swanlog.warning("You are logging a duplicate and same instance, this is not recommended")
        self.__file_name = value


class ParseResult:
    """
    数据解析结果
    """

    def __init__(
        self,
        section: str = None,
        chart: BaseType.Chart = None,
        data: Union[List[str], float] = None,
        config: Optional[List[Dict]] = None,
        more: Optional[List[Dict]] = None,
        buffers: Optional[List[MediaBuffer]] = None,
    ):
        """
        :param section: 转换后数据对应的section
        :param chart: 转换后数据对应的图表类型，枚举类型
        :param data: 存储在.log中的数据
        :param config: 存储在.log中的配置
        :param more: 存储在.log中的更多信息
        :param buffers: 存储于media文件夹中的原始数据，比特流，特别的，对于某些字符串即原始数据的情况，此处为None
        """
        self.__data = data
        self.config = config
        self.more = more
        self.buffers = buffers
        self.section = section
        self.chart = chart
        # 默认的reference
        self.reference = "step"
        self.step = None

    @property
    def data(self):
        """
        获取数据，如果不确定为何种数据类型，可以使用此属性
        """
        return self.__data

    @property
    def float(self):
        """
        当确认数据为float时，可以使用此属性明确为float
        """
        return self.__data

    @float.setter
    def float(self, value):
        if not isinstance(value, (int, float)):
            # 无法被nan或者inf表示的数据，不允许被设置为float
            if not (math.isnan(float(value)) or math.isinf(float(value))):
                raise TypeError(f"Expected float, but got {type(value)}")
            elif math.isnan(float(value)):
                swanlog.warning(f"Your are setting a nan to float")
            else:
                swanlog.warning(f"Your are setting a inf to float")
        self.__data = value

    @property
    def strings(self):
        """
        当确认数据为字符串时，可以使用此属性明确为字符串
        """
        return self.__data

    @strings.setter
    def strings(self, value):
        if not isinstance(value, list) or not all(isinstance(i, str) for i in value):
            raise TypeError(f"Expected list(str), but got {type(value)}")
        # 不允许设置为空列表，目的是希望在设置时就能保证数据的完整性
        if not len(value):
            raise ValueError("Expected list(str) with at least one element")
        self.__data = value
