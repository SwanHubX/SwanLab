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
import wave


class Audio(BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        print("root_dir", self.settings.root_dir)

        self.image = self.preprocess(self.value)
        save_path = os.path.join(self.settings.static_dir, f"image-{self.tag}-{self.step}.png")

        # 保存数据到指定目录
        self.save(save_path)
        print("save_path:", save_path)

        # 获得目录的相对路径
        save_relative_path = os.path.relpath(save_path, self.settings.root_dir)
        print("save_relative_path", save_relative_path)

        return save_relative_path

    def preprocess(self, data):
        """
        将音频文件读取为numpy array
        """
        return None

    def save(self, save_path):
        """
        保存 wav 文件到指定路径
        """
        try:
            self.data.save(save_path)
        except Exception as e:
            raise ValueError(f"Could not save the audio file to the path: {save_path}") from e

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Audio"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.line
