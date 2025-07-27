from typing import List, Union

from swanlab.toolkit import BaseType, MediaBuffer, MediaType
from .audio import Audio
from .custom_charts import echarts, Echarts, PyEchartsBase, PyEchartsTable
from .image import Image
from .line import FloatConvertible, Line
from .object3d import Object3D, Molecule
from .text import Text
from .line import Line, FloatConvertible
from typing import Union, List
from .video import Video
from .wrapper import DataWrapper
from .custom_charts.metrics import roc_curve, pr_curve, confusion_matrix

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
    "Video",
    # 模块
    "echarts",
    # metrics,
    "roc_curve",
    "pr_curve",
    "confusion_matrix",
]
