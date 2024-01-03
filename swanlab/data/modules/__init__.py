from .base import BaseType
from typing import Protocol, Union


class FloatConvertible(Protocol):
    def __float__(self) -> float:
        ...


DataType = Union[float, FloatConvertible, int, BaseType]
