"""
@author: caddiesnew
@file: html.py
@time: 2026/6/8
@description: HTML 数据类型处理模块
"""

from pathlib import Path
from typing import IO, TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab.sdk.internal.run.transforms.html import Html


HtmlDataType = Union[
    "Html",  # 套娃：嵌套 Html 对象
    str,  # 原始 HTML 字符串 或 文件路径
    Path,  # pathlib.Path 文件路径
    IO,  # 文件类对象 (TextIO)
]

HtmlDatasType = Union["Html", List["Html"], str, List[str], Path, List[Path], IO, List[IO]]
