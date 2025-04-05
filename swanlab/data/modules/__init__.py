from typing import List, Union

from swankit.core.data import BaseType, MediaBuffer, MediaType

from .audio import Audio
from .image import Image
from .line import FloatConvertible, Line
from .object3d import Model3D, Object3D, PointCloud, Molecule
from .text import Text
from .wrapper import DataWrapper

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
    "MediaBuffer",
    "Object3D",
    "PointCloud",
    "Model3D",
    "Molecule",
]
