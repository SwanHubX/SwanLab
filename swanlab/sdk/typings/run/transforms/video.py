"""
@author: cunyue
@file: video.py
@time: 2026/3/17 19:36
@description: 视频处理模块
"""

from io import BytesIO
from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab.sdk.internal.run.transforms.video import Video

VideoDataType = Union["Video", str, bytes, BytesIO]


VideoDatasType = Union["Video", List["Video"], str, List[str], bytes, List[bytes], BytesIO, List[BytesIO]]
