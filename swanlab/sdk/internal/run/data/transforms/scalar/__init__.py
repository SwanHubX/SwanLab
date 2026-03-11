"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 16:32
@description: 标量处理模块
"""

import math
from typing import Union

from swanlab.proto.swanlab.data.v1.scalar_pb2 import ScalarValue
from swanlab.sdk.internal.run.data.transforms.abc import TransformType


class Scalar(TransformType):
    def __init__(self):
        # 标量类型直接用字符串、数字、布尔值表示，不应该被实例化
        raise NotImplementedError("Scalar should not be instantiated directly.")

    @staticmethod
    def transform(data: Union[float, int, str, bool]) -> ScalarValue:
        """
        处理标量数据，包括数字、字符串、布尔值等。

        1. 布尔值：转为 1.0 或 0.0
        2. 数字（int/float）：正常转换为 float
        3. 字符串：尝试转换为数字，如果失败则使用 NaN

        :param data: 待处理的数据

        :raises TypeError: 如果数据类型不支持转换为浮点数
        """
        # 0. 鸭子类型检测：如果是 Tensor 或 numpy array，尝试提取其标量值
        if hasattr(data, "item") and callable(data.item):  # type: ignore
            try:
                # .item() 会将 0维 或 1元素张量 转化为原生 Python int/float/bool
                data = data.item()  # type: ignore
            except ValueError as e:
                # 如果传入的是多维张量 (例如 tensor([1.0, 2.0]))，.item() 会抛出 ValueError
                raise ValueError(f"无法将多元素 Tensor/Array 转换为单一标量记录。异常: {e}")

        # 1. 优先判断 bool，因为 bool 是 int 的子类
        if isinstance(data, bool):
            return ScalarValue(number=float(data))

        # 2. 判断纯数字 (显式指定 int 和 float)
        if isinstance(data, (int, float)):
            return ScalarValue(number=float(data))

        # 3. 判断字符串
        if isinstance(data, str):
            try:
                return ScalarValue(number=float(data))
            except ValueError:
                # 转换失败时，安全地返回 NaN
                return ScalarValue(number=math.nan)

        # 兜底：类型不匹配
        raise TypeError(f"Unsupported scalar type: {type(data).__name__}")
