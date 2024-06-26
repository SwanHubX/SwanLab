#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 02:34
@File: __init__.py
@IDE: pycharm
@Description:
    音频模块
"""
from swankit.core import MediaBuffer, DataSuite as D, MediaType
from typing import Union, Any

try:
    # noinspection PyPackageRequirements
    import soundfile as sf

    # noinspection PyPackageRequirements
    import numpy as np

    InputType = Union[str, np.ndarray]
except ImportError:
    sf, np = None, None
    InputType = Union[str, Any]


class Audio(MediaType):
    SF_SUPPORT_DTYPE = [np.dtype(d) for d in ["float32", "float64", "int16", "int32"]] if sf else []

    def __init__(self, data_or_path: InputType, sample_rate: int = 44100, caption: str = None):
        """Audio class constructor

        Parameters
        ----------
        data_or_path: str or numpy.ndarray
            Path to an audio file or numpy array of audio data.
        sample_rate: int
                Sample rate of the audio data. Required when input is numpy array.
        caption: str
            Caption for the audio.
        """
        if sf is None or np is None:
            raise ImportError(
                "soundfile and numpy are required for Audio class, "
                "you can install them by `pip install soundfile numpy`"
            )
        super().__init__()
        if isinstance(data_or_path, str):
            # 如果输入为路径字符串
            try:
                audio_data, sample_rate = sf.read(data_or_path)
                # 转置为 (num_channels, num_frames) 的形式
                audio_data = audio_data.T
            except Exception as e:
                raise ValueError(f"Invalid audio path: {data_or_path}") from e
        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy array ，要求输入为 (num_channels, num_frames) 的形式
            # 支持单声道 或 双声道 两种形式

            if data_or_path.dtype not in self.SF_SUPPORT_DTYPE:
                e = (
                    f"Invalid numpy array for the audio data, support dtype is {self.SF_SUPPORT_DTYPE}, "
                    f"but got {data_or_path.dtype}"
                )
                raise TypeError(e)

            # 如果data_or_path是一维, 则reshape为2维
            if len(data_or_path.shape) == 1:
                data_or_path = data_or_path.reshape(1, -1)

            # 获取通道数
            num_channels = data_or_path.shape[0]

            if num_channels != 2 and num_channels != 1:
                raise TypeError("Invalid numpy array for the audio data, support shape is (num_channels, num_frames)")
            if sample_rate is None:
                raise TypeError("sample_rate must be provided when input is numpy array while constructing Audio()")
            audio_data = data_or_path
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported audio type. Please provide a valid path or numpy array.")
        self.audio_data = audio_data
        """
        转换为矩阵后的数据
        """
        self.buffer = MediaBuffer()
        sf.write(self.buffer, audio_data.T, sample_rate, format="wav")
        self.sample_rate = sample_rate
        self.caption = D.check_caption(caption)

    # ---------------------------------- 覆写方法 ----------------------------------
    def parse(self):
        # 文件名称
        hash_name = D.get_hash_by_ndarray(self.audio_data)[:16]
        save_name = f"audio-step{self.step}-{hash_name}.wav"
        return save_name, self.buffer

    def get_more(self):
        """返回more数据"""
        # 如果传入的是Audio类列表
        return {"caption": self.caption} if self.caption is not None else None

    def get_section(self):
        """设定分组名"""
        return "Audio"

    def get_chart(self):
        """设定图表类型"""
        return self.Chart.AUDIO
