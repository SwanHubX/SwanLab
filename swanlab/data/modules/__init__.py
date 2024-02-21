from .base import BaseType
from .audio import Audio
from .image import Image
from .text import Text
from .object_3d import Object3D
from typing import Protocol, Union


class FloatConvertible(Protocol):
    def __float__(self) -> float: ...


DataType = Union[float, FloatConvertible, int, BaseType]
