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
import logging
import numpy as np
from io import BytesIO
from typing import  Optional
from moviepy.editor import VideoClip
import numpy as np


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
        data_or_path: Union[ str,"np.ndarray", "BytesIO",List["Video"]],
        caption: Optional[str] = None,
        fps: int = 4,
        format: Optional[str] = None,
        ):
        super().__init__(data_or_path)
        self.caption = self.__convert_caption(caption)
        self.fps = fps
        self.format = format
        self.time = None
        self.height = None
        self.width = None
        self.channels = None
        if self.format not in Video.EXTS:
            raise ValueError("swanlab.Video accepts %s formats" % ", ".join(Video.EXTS))

    def get_data(self, *args, **kwargs):
        # 如果传入的是Video类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        # 视频预处理
        self.video_data = self.__preprocess(self.value)
        # 获取视频hash值
        hash_name = (
            get_file_hash_numpy_array(self.video_data)[:16]
            if isinstance(self.video_data, np.ndarray)
            else get_file_hash_path(self.video_data)[:16]
        )
        # 设置视频保存路径, 保存文件名
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = (
            f"{self.caption}-step{self.step}-{hash_name}.mp4"
            if self.caption is not None
            else f"video-step{self.step}-{hash_name}.mp4"
        )
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存视频数据到指定目录
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str","np.ndarray", "BytesIO"]

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

    def __preprocess(self,data):
        """将不同类型的输入转换为"""
        if isinstance(data,str):
            # 如果输入为字符串
           self.__load_video_from_path(data)
        elif isinstance(data,np.ndarray):
             # 如果输入为numpy array
            tensor = self.__prepare_video(self.data)
            self.time, self.height, self.width, self.channels = tensor.shape
            return tensor
        elif isinstance(data,BytesIO):
            # 如果输入为二进制数据流
            self.__load_video_from_BytesIO()
            pass
        else:
            # 以上都不是，则报错
            raise ValueError(
                    "swanlab.Video accepts a file path or numpy like data as input"
                )

    def __load_video_from_path(self,path):
        """判断字符串是否为正确的视频路径，如果是则返回np.ndarry类型对象，如果不是则报错"""
        pass

    def __load_video_from_BytesIO(self,Bytes):
        """将字节流保存为临时视频文件，再读取为np.ndarry类型对象"""

        pass

    def __prepare_video(self, video: "np.ndarray") -> "np.ndarray":
        """This logic was mostly taken from tensorboardX."""
        if video.ndim < 4:
            raise ValueError(
                "Video must be atleast 4 dimensions: time, channels, height, width"
            )
        if video.ndim == 4:
            video = video.reshape(1, *video.shape)
        b, t, c, h, w = video.shape

        if video.dtype != np.uint8:
            logging.warning("Converting video data to uint8")
            video = video.astype(np.uint8)

        def is_power2(num: int) -> bool:
            return num != 0 and ((num & (num - 1)) == 0)

        # pad to nearest power of 2, all at once
        if not is_power2(video.shape[0]):
            len_addition = int(2 ** video.shape[0].bit_length() - video.shape[0])
            video = np.concatenate(
                (video, np.zeros(shape=(len_addition, t, c, h, w))), axis=0
            )

        n_rows = 2 ** ((b.bit_length() - 1) // 2)
        n_cols = video.shape[0] // n_rows

        video = np.reshape(video, newshape=(n_rows, n_cols, t, c, h, w))
        video = np.transpose(video, axes=(2, 0, 4, 1, 5, 3))
        video = np.reshape(video, newshape=(t, n_rows * h, n_cols * w, c))
        return video

    def __save(self, save_path):
        """将视频保存到指定路径"""
        fps = self.fps
        write_video_data = self.video_data
        # 计算视频的总时长
        frame_count = write_video_data.shape[0]
        duration = frame_count / fps
        # 读取视频帧
        clip = VideoClip(write_video_data, duration=duration)
        try:
            clip.write_videofile(save_path, fps=fps, codec="libx264")
        except Exception as e:
            raise TypeError(f"Could not save the video to the path: {save_path}") from e

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        # 如果传入的是Video类列表
        if isinstance(self.value, list):
            return self.get_more_list()
        else:
            return (
                {
                    "caption": self.caption,
                }
                if self.caption is not None
                else None
            )


    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Video"

    def get_chart_type(self, *args, **kwargs) -> str:
        """设定图表类型"""
        return self.chart.video

