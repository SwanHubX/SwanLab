#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-20 17:16:44
@File: swanlab\data\modules\video.py
@IDE: vscode
@Description:
     视频数据解析
"""
from .base import BaseType
from ._utils import get_file_hash_numpy_array, get_file_hash_path
import os
from typing import Union, List
import logging
import numpy as np
from io import BytesIO
from typing import Optional
import shutil


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
            The time axis is measured in frames
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
        data_or_path: Union[str, "np.ndarray", "BytesIO", List["Video"]],
        caption: Optional[str] = None,
        fps: int = 4,
        format: Optional[str] = None,
    ):
        super().__init__(data_or_path)
        self.caption = self.__convert_caption(caption)
        self.fps = fps
        self.format = format
        self.height = None
        self.width = None
        self.channels = None
        # self.format的默认值为mp4
        self.format = format or "mp4"
        # 如果self.format不在Video.EXTS中，则报错
        if self.format not in Video.EXTS:
            raise ValueError("swanlab.Video 'format' paramters's format range is ({}) .".format(", ".join(Video.EXTS)))

    def get_data(self, *args, **kwargs):
        # 如果传入的是Video类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        # 视频预处理
        self.video_data = self.__preprocess(self.value)

        # 判断是否要保存(mode='disabled'时不保存)
        if not self.settings.should_save:
            return

        # 获取视频hash值
        hash_name = (
            get_file_hash_numpy_array(self.video_data)[:16]
            if isinstance(self.video_data, np.ndarray)
            else get_file_hash_path(self.video_data)[:16]
        )
        # 设置视频保存路径, 保存文件名
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"text-step{self.step}-{hash_name}.{self.format}"
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存视频数据到指定目录
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "np.ndarray", "BytesIO"]

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
        return caption.strip() if caption else None

    def __preprocess(self, data):
        """将不同类型的输入转换为np.ndarray类型"""
        if isinstance(data, (str, BytesIO)):
            # 如果输入为字符串或二进制数据流
            return data

        elif isinstance(data, np.ndarray):
            # 如果输入为numpy array
            data = self.__numpy_preprocess(data)
            return data

        else:
            # 以上都不是，则报错
            raise ValueError("swanlab.Video accepts a file path or numpy like data as input")

    def __numpy_preprocess(self, video: "np.ndarray") -> "np.ndarray":
        """This logic was mostly taken from tensorboardX."""
        if video.ndim < 4:
            raise ValueError("Video must be atleast 4 dimensions: time, channels, height, width")
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
            video = np.concatenate((video, np.zeros(shape=(len_addition, t, c, h, w))), axis=0)

        n_rows = 2 ** ((b.bit_length() - 1) // 2)
        n_cols = video.shape[0] // n_rows

        video = np.reshape(video, newshape=(n_rows, n_cols, t, c, h, w))
        video = np.transpose(video, axes=(2, 0, 4, 1, 5, 3))
        video = np.reshape(video, newshape=(t, n_rows * h, n_cols * w, c))

        _, self.height, self.width, self.channels = video.shape
        return video

    def __save(self, save_path):
        """将视频保存到指定路径"""
        # 如果data_or_path是路径, 则直接复制文件到目标文件夹
        if isinstance(self.video_data, str):
            shutil.copyfile(self.video_data, save_path)

        # 如果data_or_path是字节流, 则直接写入文件到目标文件夹
        elif isinstance(self.video_data, BytesIO):
            with open(save_path, "wb") as f:
                f.write(self.video_data.read())

        # 如果data_or_path是np.ndarray，则使用moviepy库保存视频
        elif isinstance(self.video_data, np.ndarray):
            try:
                import moviepy
            except ImportError as e:
                raise TypeError(
                    "swanlab.Video requires moviepy when process numpy data. Install with 'pip install moviepy'."
                )

            from moviepy.editor import ImageSequenceClip

            write_video_data_bgr = self.video_data
            # OpenCV 默认读取为BGR，转化为RGB
            write_video_data_rgb = [frame[:, :, ::-1] for frame in write_video_data_bgr]
            # 读取视频帧
            clip = ImageSequenceClip(list(write_video_data_rgb), fps=self.fps)
            try:
                clip.write_videofile(save_path, codec="libx264")
                clip.close()
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
