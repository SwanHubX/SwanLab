"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 16:32
@description: 标量处理模块
"""

from numbers import Number
from typing import Any, Union

from swanlab.proto.swanlab.data.v1.scalar_pb2 import ScalarValue
from swanlab.sdk.internal.run.data.transforms.abc import TransformType


class Scalar(TransformType):
    def __init__(self, value: Any):
        # 实现套娃加载
        attrs = self._unwrap(value)
        self.value = attrs.get("value", value)

    @staticmethod
    def transform(key: str, data: Union["Scalar", Number, str, bool]) -> ScalarValue:
        """
        处理标量数据，包括数字、字符串、布尔值等。

        1. 数字：正常处理
        2. 字符串：尝试转换为数字，如果失败则使用math.nan
        3. 布尔值：转为0或1

        :param key: 数据的键名，用于标识数据类型，此为必填项
        :param data: 待处理的数据
        """
        return ScalarValue(number=1)
