"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 13:06
@description: SwanLab 数据转换模块，将用户输入的数据封装为Protobuf格式
"""

from typing import Any, List, Type, Union

from swanlab.sdk.internal.context import TransformMedia

from .audio import Audio
from .image import Image
from .scalar import Scalar
from .text import Text

__all__ = ["Text", "Scalar", "Audio", "Image", "normalize_media_input"]


def normalize_media_input(
    media_cls: Type[TransformMedia],
    data: Union[Any, List[Any]],
    **kwargs,
) -> List[TransformMedia]:
    """
    规范化媒体输入为统一的列表格式。

    :param media_cls: 媒体类型类（如 Text, Image 等）
    :param data: 单个数据或数据列表（可以是原始数据或已实例化的媒体对象）
    :param kwargs: 额外参数（如 caption），可以是单个值或列表
    :return: 媒体对象列表
    """
    # 如果 data 已经是目标类型的实例列表且无额外参数，直接返回
    if isinstance(data, list) and data and isinstance(data[0], media_cls) and not kwargs:
        return data

    # 如果 data 是单个目标类型实例且无额外参数，返回列表
    if isinstance(data, media_cls) and not kwargs:
        return [data]

    # 规范化为列表
    data_list = data if isinstance(data, list) else [data]
    n = len(data_list)

    # 规范化 kwargs：单个值扩展为列表，列表保持不变
    normalized_kwargs = {}
    for key, value in kwargs.items():
        if value is None:
            normalized_kwargs[key] = [None] * n
        elif isinstance(value, list):
            if len(value) != n:
                raise ValueError(f"Length of {key} ({len(value)}) must match data length ({n})")
            normalized_kwargs[key] = value
        else:
            normalized_kwargs[key] = [value] * n

    # 构造媒体对象列表
    result = []
    for i, item in enumerate(data_list):
        if isinstance(item, media_cls):
            result.append(item)
        else:
            item_kwargs = {k: v[i] for k, v in normalized_kwargs.items()}
            result.append(media_cls(*[item], **item_kwargs))

    return result
