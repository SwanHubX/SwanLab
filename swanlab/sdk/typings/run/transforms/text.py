"""
@author: cunyue
@file: text.py
@time: 2026/3/17 19:36
@description: 文本处理模块
"""

from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab.sdk.internal.run.transforms.text import Text


TextDataType = Union[
    "Text",
    str,
]

TextDatasType = Union["Text", List["Text"], str, List[str]]
