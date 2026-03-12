"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 16:32
@description: 标量处理模块
"""

import math
from typing import Any, Union

from swanlab.proto.swanlab.data.v1.scalar_pb2 import ScalarValue
from swanlab.sdk.internal.context import TransformType
from swanlab.sdk.utils.helper import catch_and_return_none


class Scalar(TransformType):
    def __init__(self):
        # 标量类型直接用字符串、数字、布尔值表示，不应该被实例化
        raise NotImplementedError("Scalar should not be instantiated directly.")

    @staticmethod
    def transform(data: Any) -> ScalarValue:
        """
        处理标量数据，包括数字、字符串、布尔值等。

        1. 布尔值：转为 1.0 或 0.0
        2. 数字（int/float）：正常转换为 float
        3. 字符串：尝试转换为数字，如果能转换为数字则转换为数字，否则报错，额外处理 NaN 和 Inf 的情况，此时都返回 NaN

        :param data: 待处理的数据

        :raises TypeError: 如果数据类型不支持转换为浮点数
        """
        # 0. 鸭子类型检测：如果是 Tensor 或 numpy array，尝试提取其标量值
        this_value = _transform_tensor_or_array(data)
        if this_value is None:
            full_type_name = f"{type(data).__module__}.{type(data).__name__}"
            raise TypeError(
                f"Failed to extract scalar value from {full_type_name}. "
                "If it's a Tensor or Array, please ensure it's a scalar value.",
            )
        data = this_value
        # 1. 优先判断 bool，因为 bool 是 int 的子类
        if isinstance(data, bool):
            return ScalarValue(number=float(data))

        # 2. 判断纯数字 (显式指定 int 和 float)
        if isinstance(data, (int, float)):
            return ScalarValue(number=float(data))

        # 3. 判断字符串
        if isinstance(data, str):
            try:
                value = float(data)
            except ValueError:
                raise TypeError(f"Unsupported scalar string value: '{data}'.")
            # NaN 和 Inf 都归一化为 NaN
            if math.isnan(value) or math.isinf(value):
                return ScalarValue(number=math.nan)
            return ScalarValue(number=value)

        # 兜底：类型不匹配
        raise TypeError(f"Unsupported scalar type: {type(data).__name__}.")


@catch_and_return_none()
def _transform_tensor_or_array(data: Any) -> Union[float, int, str, bool]:
    if hasattr(data, "item") and callable(data.item):
        return data.item()  # type: ignore
    return data
