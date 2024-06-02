# -*- coding: utf-8 -*-
"""
Date: 2024-01-22 15:01:36
IDE: VSCode
FilePath: /SwanLab/swanlab/data/modules/audio.py
Description:
    音频数据解析
"""
from .base import BaseType
from ._utils import get_file_hash_numpy_array, get_file_hash_path
import os
from typing import Union, List

### 以下为音频数据解析的依赖库
import soundfile as sf
import numpy as np
import json
import os


class Audio(BaseType):
    """Audio class constructor

    Parameters
    ----------
    data_or_path: str or numpy.array
        Path to an audio file or numpy array of audio data.
    sample_rate: int
            Sample rate of the audio data. Required when input is numpy array.
    caption: str
        Caption for the audio.
    """

    def __init__(
        self,
        data_or_path: Union[str, np.ndarray, List["Audio"]],
        sample_rate: int = 44100,
        caption: str = None,
    ):
        """Accept a path to an audio file on a numpу array of audio data."""
        super().__init__(data_or_path)
        self.audio_data = None
        self.sample_rate = sample_rate
        self.caption = self.__convert_caption(caption)

    def get_data(self):
        # 如果传入的是Audio类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        self.__preprocess(self.value)
        # 判断是否要保存(mode='disabled'时不保存)
        if not self.settings.should_save:
            return
        hash_name = (
            get_file_hash_numpy_array(self.audio_data)[:16]
            if isinstance(self.audio_data, np.ndarray)
            else get_file_hash_path(self.audio_data)[:16]
        )
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"audio-step{self.step}-{hash_name}.wav"
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存音频数据到指定目录
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "numpy.ndarray"]

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

    def __preprocess(self, data_or_path):
        """
        根据输入不同的输入类型进行不同处理
        """
        if isinstance(data_or_path, str):
            # 如果输入为路径字符串
            self.audio_data, self.sample_rate = self.__load_audio_from_path(data_or_path)

        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy array ，要求输入为 (num_channels, num_frames) 的形式
            # 支持单声道 或 双声道 两种形式

            SF_SUPPORT_DTYPE = [np.dtype(d) for d in ["float32", "float64", "int16", "int32"]]

            if data_or_path.dtype not in SF_SUPPORT_DTYPE:
                raise TypeError(
                    f"Invalid numpy array for the audio data, support dtype is {SF_SUPPORT_DTYPE}, but got {data_or_path.dtype}"
                )

            # 如果data_or_path是一维, 则reshape为2维
            if len(data_or_path.shape) == 1:
                data_or_path = data_or_path.reshape(1, -1)

            # 获取通道数
            num_channels = data_or_path.shape[0]

            if num_channels != 2 and num_channels != 1:
                raise TypeError("Invalid numpy array for the audio data, support shape is (num_channels, num_frames)")
            if self.sample_rate is None:
                raise TypeError("sample_rate must be provided when input is numpy array while constructing Audio()")
            self.audio_data = data_or_path
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported audio type. Please provide a valid path or numpy array.")

    def __load_audio_from_path(self, path):
        """判断字符串是否为正确的音频文件路径，如果是则转换为numpy array，如果不是则报错"""
        try:
            audio_data, sample_rate = sf.read(path)
            # 转置为 (num_channels, num_frames) 的形式
            audio_data = audio_data.T
            return audio_data, sample_rate

        except Exception as e:
            raise ValueError(f"Invalid audio path: {path}") from e

    def __save(self, save_path):
        """
        保存媒体资源文件 .wav 到指定路径
        audio-step{}.wav
        """

        try:
            write_audio_data = self.audio_data.T  # type: np.ndarray
            audio_bytes = write_audio_data.tobytes()  # 字节流

            sf.write(save_path, write_audio_data, self.sample_rate)
        except Exception as e:
            raise ValueError(f"Could not save the audio file to the path: {save_path}") from e

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        # 如果传入的是Audio类列表
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
        return "Audio"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.audio
