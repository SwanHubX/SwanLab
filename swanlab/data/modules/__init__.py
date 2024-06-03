from .base import BaseType, MediaType, ParseResult
from .audio import Audio
from .image import Image
from .text import Text
from .line import Line, FloatConvertible
from typing import Union, List
from .wrapper import DataWrapper, ErrorInfo

DataType = Union[int, float, FloatConvertible, BaseType, List[BaseType]]
ChartType = BaseType.Chart

__all__ = [
    "FloatConvertible",
    "DataWrapper",
    "MediaType",
    "Image",
    "Audio",
    "Text",
    "Line",
    "DataType",
    "ChartType",
    "ErrorInfo"
]
