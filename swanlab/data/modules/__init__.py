from .base import BaseType
from .audio import Audio
from .image import Image
from .text import Text
from typing import Protocol, Union, List


class FloatConvertible(Protocol):
    def __float__(self) -> float: ...


DataType = Union[float, FloatConvertible, int, BaseType, List[BaseType]]
