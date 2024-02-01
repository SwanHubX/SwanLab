# -*- coding: utf-8 -*-
"""
Date: 2024-01-22 15:01:36
IDE: VSCode
FilePath: /SwanLab/swanlab/data/modules/audio.py
Description:
    音频数据解析
"""
from .base import BaseType
import os

### 以下为音频数据解析的依赖库
import soundfile as sf
import numpy as np

# import matplotlib.pyplot as plt
from IPython.display import Audio as AudioDisplay
import json


class Audio(BaseType):
    """
    Audio class constructor
    ----
    data_or_path: str or numpy.array
    """

    def __init__(self, data_or_path, sample_rate: int = None):
        """Accept a path to an audio file on a numpу array of audio data."""
        super().__init__(data_or_path)
        self.audio_data = None
        self.sample_rate = None
        if sample_rate is not None:
            self.sample_rate = sample_rate

    def get_data(self):
        # print(self.step, self.tag, self.settings.static_dir)

        # numpy.array() 类型保存音频数据
        self.preprocess(self.value)
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"audio-{self.step}.wav"
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存音频数据到指定目录
        self.save(save_path)
        # print("save_path:", save_path)
        # print("save_relative_path", save_relative_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "numpy.array"]

    def preprocess(self, data_or_path):
        """
        根据输入不同的输入类型进行不同处理
        """
        # print("data_or_path", data_or_path)
        if isinstance(data_or_path, str):
            # 如果输入为路径字符串
            # 根据输入是否为 json , 选择不同的加载方式
            if data_or_path.endswith(".json"):
                self.load_from_json(data_or_path)
            else:
                self.audio_data, self.sample_rate = self.load_audio_from_path(data_or_path)

        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy array ，要求输入为 (num_channels, num_frames) 的形式
            # 支持单声道 或 双声道 两种形式
            num_channels = data_or_path.shape[0]
            if num_channels != 2 and num_channels != 1:
                raise ValueError("Invalid numpy array for the audio data, support shape is (num_channels, num_frames)")
            if self.sample_rate is None:
                raise ValueError("sample_rate must be provided when input is numpy array while constructing Audio()")
            self.audio_data = data_or_path
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported audio type. Please provide a valid path or numpy array.")

    def load_audio_from_path(self, path):
        """判断字符串是否为正确的音频文件路径，如果是则转换为numpy array，如果不是则报错"""
        try:
            audio_data, sample_rate = sf.read(path)
            # 转置为 (num_channels, num_frames) 的形式
            audio_data = audio_data.T
            return audio_data, sample_rate

        except Exception as e:
            raise ValueError(f"Invalid audio path: {path}") from e

    # def plot_spectrogram(self, plot_save_path):
    #     """获取音频的频谱"""
    #     num_channels = self.audio_data.shape[0]
    #     sample_rate = self.sample_rate
    #     figure, axes = plt.subplots(num_channels, 1)
    #     if num_channels == 1:
    #         axes = [axes]
    #     for c in range(num_channels):
    #         axes[c].specgram(self.audio_data[c], Fs=sample_rate)
    #         if num_channels > 1:
    #             axes[c].set_ylabel(f"Channel {c+1}")
    #     figure.savefig(plot_save_path)
    #     plt.close()

    def to_json(self, audio_path, json_path):
        """将音频元数据转换为json文件"""
        audio_json = {}
        audio_json["type"] = "audio"
        audio_json["file_path"] = audio_path
        audio_json["sample_rate"] = self.sample_rate
        audio_json["duration"] = self.audio_data.shape[1]

        with open(json_path, "w") as f:
            json.dump(audio_json, f)

    def load_from_json(self, json_path):
        """从json文件中加载音频元数据"""
        with open(json_path, "r") as f:
            audio_json = json.load(f)
        json_sample_rate = audio_json["sample_rate"]
        file_path = audio_json["file_path"]
        audio_data, data_sample_rate = self.load_audio_from_path(file_path)
        if json_sample_rate != data_sample_rate:
            raise ValueError("sample_rate in json file is not equal to the sample_rate of the audio file")
        self.audio_data = audio_data
        self.sample_rate = data_sample_rate

    def get_play_audio(self):
        """获取可播放音频类型"""
        try:
            return AudioDisplay(self.audio_data, rate=self.sample_rate)
        except Exception as e:
            raise ValueError(f"Could not play the audio file") from e

    def save(self, save_path):
        """
        保存静态资源文件 .wav 到指定路径
        audio-{tag}-{step}.wav
        audio-{tag}-{step}.json
        """
        try:
            write_audio_data = self.audio_data.T
            sf.write(save_path, write_audio_data, self.sample_rate, subtype="PCM_16")

            # 这边存 json 的时候放的是相对路径
            # json_path = save_path.replace(".wav", ".json")
            # save_relative_path = os.path.relpath(save_path, self.settings.root_dir)
            # self.to_json(save_relative_path, json_path)

            # 保存频谱图
            # spectrogram_save_path = save_path.replace(".wav", "-spectrogram.png")
            # self.plot_spectrogram(spectrogram_save_path)

        except Exception as e:
            raise ValueError(f"Could not save the audio file to the path: {save_path}") from e

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Audio"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.audio
