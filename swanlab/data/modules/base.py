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
from ..settings import SwanDataSettings
from typing import List, Dict, Optional, ByteString, Union, Tuple, Callable
from enum import Enum
import hashlib
import math
import io


class U:
    """
    BaseType在处理时的通用工具类
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
        # 当前运行时配置
        self.settings: Optional[SwanDataSettings] = None

    def inject(self, key: str, step: int, settings: SwanDataSettings):
        """
        注入属性
        :param key: key名称
        :param step: 当前步骤
        :param settings: 当前运行时配置
        """
        self.key = key
        self.step = step
        self.settings = settings


class BaseType(ABC, DynamicProperty, U):
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

            def __init__(self, chart_type: str):
                self.chart_type = chart_type

        LINE = ChartItem("line")

        IMAGE = ChartItem("image")

        AUDIO = ChartItem("audio")

        TEXT = ChartItem("text")

    # ---------------------------------- 需要子类实现的方法 ----------------------------------

    @abstractmethod
    def parse(self) -> Tuple[Union[str, float], Optional[ByteString]]:
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

    def get_raw_write_handler(self) -> Optional[Callable[[str], None]]:
        """
        获取写入raw的处理函数，约定处理函数的输入仅有一个参数，即存储路径
        处理函数内不需要对路径做处理，调用时会处理
        没有特定的写入函数则不需要覆写
        """
        pass


class MediaType(BaseType):  # noqa
    """
    媒体类型，用于区分标量和媒体，不用做实例化，应该由子类继承
    """
    pass


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
        raw: Optional[List[ByteString]] = None,
        raw_write_handler: Optional[Callable[[str], None]] = None,
    ):
        """
        :param section: 转换后数据对应的section
        :param chart: 转换后数据对应的图表类型，枚举类型
        :param data: 存储在.log中的数据
        :param config: 存储在.log中的配置
        :param more: 存储在.log中的更多信息
        :param raw: 存储于media文件夹中的原始数据，比特流，特别的，对于某些字符串即原始数据的情况，此处为None
        :param raw_write_handler: 写入raw的处理函数
        """
        self.data = data
        self.config = config
        self.more = more
        self.raw = raw
        self.section = section
        self.chart = chart
        self.write_handler = raw_write_handler
        # 默认的reference
        self.reference = "step"
