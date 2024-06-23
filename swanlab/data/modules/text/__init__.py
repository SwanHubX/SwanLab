#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 02:35
@File: __init__.py
@IDE: pycharm
@Description:
    文本模块
"""
from swankit.core.data import MediaType
from swankit.core import DataSuite as D
from typing import Union


class Text(MediaType):
    """Text class constructor

    Parameters
    ----------
    data: str, float, int
        text data.
    caption: str
        caption for the data, it will be displayed in the dashboard.
        e.g. swanlab.Text("Hello World", caption="This is a caption for the text data.")
    """

    def __init__(self, data: Union[str, int, float], caption: str = None):
        super().__init__()
        # 处理文本数据

        if isinstance(data, str):
            self.text_data = data
        elif isinstance(data, (int, float)):
            self.text_data = str(data)
        else:
            raise TypeError("data must be a string, int or float.")

        # 处理caption
        self.caption = D.check_caption(caption)

    # ---------------------------------- 覆写 ----------------------------------

    def parse(self):
        return self.text_data, None

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        return {"caption": self.caption} if self.caption is not None else None

    def get_section(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Text"

    def get_chart(self):
        """设定图表类型"""
        return self.Chart.TEXT
