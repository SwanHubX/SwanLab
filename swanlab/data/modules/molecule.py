# -*- coding: utf-8 -*-
"""
Date: 2024-01-22 15:01:36
IDE: VSCode
FilePath: /SwanLab/swanlab/data/modules/audio.py
Description:
    音频数据解析
"""
from .base import BaseType
from typing import Union, List


class Molecule(BaseType):
    """Molecule class constructor

    Parameters
    ----------
    data_or_path: str
    """

    def __init__(
        self,
        data_or_path: Union[str, List["Molecule"]],
        caption: str = None,
    ):

        super().__init__(data_or_path)
        self.molecule_data = None
        self.caption = self.__convert_caption(caption)

    def get_data(self):
        pass

    def expect_types(self, *args, **kwargs) -> list:
        return ["str"]

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
        pass

    def __save(self, save_path):
        pass

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
        return "Molecule"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.melocule
