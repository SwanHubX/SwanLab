"""
@author: cunyue
@file: transformer.py
@time: 2026/3/12 00:49
@description: SwanLab 运行时数据转换器抽象类，负责将用户输入的数据转换为Protobuf格式
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, List

from google.protobuf.message import Message
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord


class TransformData(ABC):
    """
    SwanLab 数据转换模块抽象基类，定义了将数据转换为Protobuf格式的方法。

    我们约定“TransformType”的相关实现子类实例仅承担“数据传输”职责，数据转换交由每个“TransformType”的实现子类的静态方法“transform”完成。

    设计为静态方法的另一个原因是考虑到大规模数据传输时，静态方法可以避免实例化开销，提高性能。

    每个子类在实现 __init__ 时必须考虑套娃问题，我们约定外层参数的优先级高于内层参数，例如：
    >>> from swanlab.proto.swanlab.metric.data.v1.media.text_pb2 import TextItem, TextValue
    >>>
    >>> class MyTransform(TransformData):
    >>>     def __init__(self, text: str | MyTransform, foo: Any = None):
    >>>         super().__init__()
    >>>         attrs = self._unwrap(text)
    >>>         self.text = attrs.get("text", text)
    >>>         self.foo = foo if foo is not None else attrs.get("foo")
    >>>     def transform(self, *, key: str, step: int, data: str = None) -> TextItem:
    >>>         return TextItem(filename=f"{key}-{step:03d}.__swanlab__.txt", content=data or "")
    >>>     @classmethod
    >>>     def build_data_record(cls,*, key: str, step: int, timestamp: Timestamp, data: Any) -> DataRecord:
    >>>         value = TextValue(items=[data])
    >>>         return DataRecord(key=key, step=step, texts=value)
    >>>     @classmethod
    >>>     def column_type(cls) -> ColumnType:
    >>>         return ColumnType.COLUMN_TYPE_TEXT

    那么此时：

    >>> # 如果内层已经定义了 foo="old"，但外层显式给了 foo="new"
    >>> t = MyTransform(MyTransform("hello", foo="old"), foo="new")
    >>> # 结果应该是 t.foo == "new"

    >>> # 生成 MetricRecord
    >>> t = MyTransform("hello")
    >>> metric_record = t.build_data_record(t.transform("key", 1))
    >>> print(metric_record)
    """

    def __init__(self, *args: Any, **kwargs: Any): ...

    @classmethod
    def _unwrap(cls, instance_or_val: Any) -> dict:
        """通用的解包辅助函数，提取实例属性，用于实现套娃加载"""
        if isinstance(instance_or_val, TransformData):
            # 返回实例的所有属性（注意排除私有属性）
            return {k: v for k, v in vars(instance_or_val).items() if not k.startswith("_")}
        return {}

    @abstractmethod
    def transform(self, *args: Any, **kwargs: Any) -> Message:
        """
        将数据转换为 Protobuf 格式。
        """
        # 套娃加载需求详见 https://github.com/SwanHubX/SwanLab/issues/1367
        ...

    @classmethod
    @abstractmethod
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: Any) -> DataRecord:
        """构建 DataRecord envelope"""
        ...

    @classmethod
    @abstractmethod
    def column_type(cls) -> ColumnType:
        """
        返回当前 TransformType 的列数据类型
        """
        ...


class TransformMedia(TransformData, ABC):
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

    @classmethod
    @abstractmethod
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[Any]) -> DataRecord:
        """构建 DataRecord envelope"""
        ...

    @abstractmethod
    def transform(self, *, step: int, path: Path) -> Message:
        """
        将媒体数据转换为Protobuf格式，并将结果写入指定目录下
        :param step: 步数
        :param path: 存储目录
        :return: Protobuf消息
        """
        ...
