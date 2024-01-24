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
import matplotlib.pyplot as plt
from IPython.display import Audio as AudioDisplay


class Audio(BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        print("root_dir", self.settings.root_dir)
        # numpy.array() 类型保存
        self.audio_data, self.sample_rate = self.preprocess(self.value)

        if os.path.exists(self.settings.static_dir) is False:
            os.makedirs(self.settings.static_dir)
        # 保存音频文件、频谱图、波形图到指定目录
        save_path = os.path.join(self.settings.static_dir, f"audio-{self.tag}-{self.step}.wav")
        self.save(save_path)
        # print("save_path:", save_path)
        # 获得目录的相对路径
        save_relative_path = os.path.relpath(save_path, self.settings.root_dir)
        print("save_relative_path", save_relative_path)
        return save_relative_path

    def preprocess(self, data_or_path, sample_rate=None):
        """
        根据输入不同的输入类型进行不同处理
        """
        if isinstance(data_or_path, str):
            # 如果输入为路径字符串
            audio_data, sample_rate = self.load_audio_from_path(data_or_path)
        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy array ，要求输入为 (num_channels, num_frames) 的形式
            # 支持单声道 或 双声道 两种形式
            if data_or_path.shape[0] != 2 or data_or_path.shape[0] != 1 or data_or_path.shape[1] < 1:
                raise ValueError("Invalid numpy array for the audio data, support shape is (num_channels, num_frames)")
            if sample_rate is None:
                raise ValueError("sample_rate must be provided when input is numpy array")
            audio_data = data_or_path
            sample_rate = sample_rate
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported audio type. Please provide a valid path or numpy array.")
        return audio_data, sample_rate

    def load_audio_from_path(self, path):
        """判断字符串是否为正确的音频文件路径，如果是则转换为numpy array，如果不是则报错"""
        try:
            audio_data, sample_rate = sf.read(path)
            # 转置为 (num_channels, num_frames) 的形式
            audio_data = audio_data.T
            return audio_data, sample_rate
        except Exception as e:
            raise ValueError(f"Invalid audio path: {path}") from e

    def plot_waveform(self, plot_save_path):
        """获取音频的波形"""
        num_channels = self.audio_data.shape[0]
        num_frames = self.audio_data.shape[1]
        sample_rate = self.sample_rate
        time_axis = np.arange(0, num_frames) / sample_rate
        figure, axes = plt.subplots(num_channels, 1)
        if num_channels == 1:
            axes = [axes]
        for c in range(num_channels):
            axes[c].plot(time_axis, self.audio_data[c], linewidth=1)
            axes[c].grid(True)
            if num_channels > 1:
                axes[c].set_ylabel(f"Channel {c+1}")
        figure.savefig(plot_save_path)
        plt.close()

    def plot_spectrogram(self, plot_save_path):
        """获取音频的频谱"""
        num_channels = self.audio_data.shape[0]
        sample_rate = self.sample_rate
        figure, axes = plt.subplots(num_channels, 1)
        if num_channels == 1:
            axes = [axes]
        for c in range(num_channels):
            axes[c].specgram(self.audio_data[c], Fs=sample_rate)
            if num_channels > 1:
                axes[c].set_ylabel(f"Channel {c+1}")
        figure.savefig(plot_save_path)
        plt.close()

    def get_play_audio(self):
        """获取可播放音频类型"""
        try:
            return AudioDisplay(self.audio_data, rate=self.sample_rate)
        except Exception as e:
            raise ValueError(f"Could not play the audio file") from e

    def save(self, save_path):
        """
        保存 wav 文件到指定路径
        绘制频谱图和波形图，保存到同样的路径
        audio-{tag}-{step}.wav
        audio-{tag}-{step}-spectrogram.png
        audio-{tag}-{step}-waveform.png
        """
        try:
            # 保存音频
            sf.write(save_path, self.audio_data, self.sample_rate, subtype="PCM_16")
            # 保存频谱图
            spectrogram_save_path = save_path.replace(".wav", "-spectrogram.png")
            self.plot_spectrogram(spectrogram_save_path)
            # 保存波形图
            waveform_save_path = save_path.replace(".wav", "-waveform.png")
            self.plot_waveform(waveform_save_path)
        except Exception as e:
            raise ValueError(f"Could not save the audio file to the path: {save_path}") from e

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Audio"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        # TODO: audio chart type 在 Chart object 中还未定义
        # 目前计划缓存 音频 + 频谱图 + 波形图
        # 可视化时，可以选择播放音频，或者显示频谱图、波形图
        return self.chart.image
