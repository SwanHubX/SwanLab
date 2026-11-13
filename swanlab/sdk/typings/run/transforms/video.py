"""
@author: cunyue
@file: video.py
@time: 2026/3/17 19:36
@description: 视频处理模块
"""

from io import BytesIO
from os import PathLike
from typing import TYPE_CHECKING, List, Literal, Protocol, Union

if TYPE_CHECKING:
    import numpy as np

    from swanlab.sdk.internal.run.transforms.video import Video

VideoFormat = Literal["gif", "mp4", "webm", "ogg"]


class NumpyConvertible(Protocol):
    """An object whose numpy() method returns uint8 TCHW or BTCHW frames."""

    def numpy(self) -> "np.ndarray": ...


VideoInput = Union["Video", str, PathLike[str], bytes, BytesIO, "np.ndarray", NumpyConvertible]
VideoDataType = VideoInput
VideoDatasType = Union[VideoInput, List[VideoInput]]
