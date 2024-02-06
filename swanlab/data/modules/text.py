# -*- coding: utf-8 -*-
from .base import BaseType
import os
from typing import Union, List


class Text(BaseType):
    """Text class constructor

    Parameters
    ----------
    data: str, float, int
        text data.
    caption: str
        caption for the data, it will be displayed in the dashboard.
        e.g. swanlab.Text("Hello World", caption="This is a caption for the text data.")
    """

    def __init__(self, data: Union[str, List["Text"]], caption: str = None):
        super().__init__(data)
        self.text_data = None
        self.caption = self.__convert_caption(caption)

    def get_data(self):
        # 如果传入的是Text类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        else:
            # 预处理文本数据
            self.__preprocess(self.value)

            # 设置文本数据的保存路径
            save_dir = os.path.join(self.settings.static_dir, self.tag)
            save_name = (
                f"{self.caption}-step{self.step}-{self.text_data[:16]}.txt"
                if self.caption is not None
                else f"text-step{self.step}-{self.text_data[:16]}.txt"
            )
            # 如果路径不存在，则创建路径
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)
            save_path = os.path.join(save_dir, save_name)

            # 保存文本数据写入到指定目录的指定json文件
            self.__save(save_path)
            return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "int", "float"]

    def __preprocess(self, data):
        """
        根据不同的输入类型进行不同处理
        """
        if isinstance(data, str):
            self.text_data = data
        elif isinstance(data, (int, float)):
            self.text_data = str(data)
        else:
            raise TypeError("data must be a string, int or float.")

    def __save(self, save_path):
        """
        将文本数据写入到指定目录的txt文件
        """
        try:
            with open(save_path, "w+", encoding="utf-8") as f:
                # 将文本数据写入到指定目录的txt文件
                f.write(self.text_data)
        except Exception as e:
            raise ValueError(f"Could not save the text data to the path: {save_path}") from e

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

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        # 如果传入的是Text类列表
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
        return "Text"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.text
