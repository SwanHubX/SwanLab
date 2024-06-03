#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 02:35
@File: __init__.py
@IDE: pycharm
@Description:
    文本模块
"""
from ..base import BaseType
from typing import Union, Tuple, Optional, ByteString


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
        self.caption = caption.strip() if caption else None

    # ---------------------------------- 覆写 ----------------------------------

    def parse(self) -> Tuple[Union[str, float], Optional[ByteString]]:
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
