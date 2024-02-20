#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-20 17:16:44
@File: swanlab\data\modules\video.py
@IDE: vscode
@Description:
    Swanlab的video数据类
"""
from .base import BaseType
from ..utils.file import get_file_hash_numpy_array, get_file_hash_path
import os
from typing import Union, List
import numpy as np
from io import BytesIO
from typing import TextIO, Optional

class Video(BaseType):
    """Video class constructor

    Parameters
    ----------
    data_or_path: numpy array, string, io
            Video can be initialized with a path to a file or an io object.
            The format must be "gif", "mp4", "webm" or "ogg".
            The format must be specified with the format argument.
            Video can be initialized with a numpy tensor.
            The numpy tensor must be either 4 dimensional or 5 dimensional.
            Channels should be (time, channel, height, width) or
            (batch, time, channel, height width)
    caption: str
        caption associated with the video for display
    fps: int
        frames per second for video. Default is 4.
    format: str
        format of video, necessary if initializing with path or io object.


    """

    EXTS = ("gif", "mp4", "webm", "ogg")
    width: Optional[int]
    height: Optional[int]


    def __init__(
        self,
        data_or_path: Union[ str,"np.ndarray", "TextIO", "BytesIO",List["Video"]],
        caption: Optional[str] = None,
        fps: int = 4,
        format: Optional[str] = None,
        ):
        super().__init__(data_or_path)
        self.caption = self.__convert_caption(caption)
        self.fps = fps
        self.format = format
        self.width = None
        self.height = None
        self.channels = None
        if self._format not in Video.EXTS:
            raise ValueError("swanlab.Video accepts %s formats" % ", ".join(Video.EXTS))

    def get_data(self, *args, **kwargs):
        return super().get_data(*args, **kwargs)





    def expect_types(self, *args, **kwargs) -> list:
        return ["str","np.ndarray", "TextIO", "BytesIO"]

    def __convert_caption(self, caption):
        """将caption转换为字符串"""
        # 如果类型是字符串，则不做转换
        if isinstance(caption, str):
            caption = caption
        # 如果类型是数字，则转换为字符串
        elif isinstance(caption, (int, float)):
            caption = str(caption)
        # 如果类型是None，则转换为默认字符串
        elif caption is None:
            caption = None
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Video"

    def get_chart_type(self, *args, **kwargs) -> str:
        """设定图表类型"""
        return self.chart.video

