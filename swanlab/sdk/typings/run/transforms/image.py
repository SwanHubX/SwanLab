"""
@author: cunyue
@file: image.py
@time: 2026/3/17 19:36
@description: 图像处理模块
"""

from typing import TYPE_CHECKING, List, Literal, Optional, Union

if TYPE_CHECKING:
    from swanlab import vendor
    from swanlab.sdk.internal.run.transforms.image import Image

ImageDataType = Union[
    "Image",
    str,
    "vendor.PIL.Image.Image",
    "vendor.np.ndarray",
    "vendor.torch.Tensor",
    "vendor.matplotlib.figure.Figure",
]

ImageDatasType = Union[
    "Image",
    List["Image"],
    str,
    List[str],
    "vendor.PIL.Image.Image",
    List["vendor.PIL.Image.Image"],
    "vendor.np.ndarray",
    List["vendor.np.ndarray"],
    "vendor.torch.Tensor",
    List["vendor.torch.Tensor"],
]

ImageModeType = Optional[str]

ImageModesType = Optional[Union[str, List[str]]]

_Files = Literal["png", "jpg", "jpeg", "bmp"]

ImageFileType = Optional[_Files]

ImageFilesType = Optional[Union[_Files, List[_Files]]]

ImageSizeType = Optional[Union[int, list, tuple]]

ImageSizesType = Optional[Union[int, List[int], tuple, List[tuple], list, List[list]]]
