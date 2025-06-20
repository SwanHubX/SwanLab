from typing import List, Union

from swanlab.toolkit import BaseType, MediaBuffer, MediaType
from .audio import Audio
from .custom_charts import echarts, Echarts, PyEchartsBase, PyEchartsTable
from .image import Image
from .line import FloatConvertible, Line
from .object3d import Object3D, Molecule
from .text import Text
from .wrapper import DataWrapper

DataType = Union[
    int,
    float,
    FloatConvertible,
    BaseType,
    List[BaseType],
    PyEchartsBase,
    PyEchartsTable,
    List[PyEchartsTable],
    List[PyEchartsBase],
]

ChartType = BaseType.Chart

__all__ = [
    # 数据类型
    "FloatConvertible",
    "DataWrapper",
    "MediaType",
    "PyEchartsBase",
    "PyEchartsTable",
    "DataType",
    "ChartType",
    "MediaBuffer",
    # 支持的图表类
    "Image",
    "Audio",
    "Text",
    "Line",
    "Object3D",
    "Molecule",
    "Echarts",
    # 模块
    "echarts",
]
