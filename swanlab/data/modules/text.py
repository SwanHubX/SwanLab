# -*- coding: utf-8 -*-
from .base import BaseType
import os
import ujson


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

    text_datas = {}

    def __init__(self, data, caption: str = None):
        super().__init__(data)
        self.text_data = None
        self.caption = caption
        if self.caption is not None:
            self.caption = self.__convert_caption(caption)

    def get_data(self):
        # 预处理文本数据
        self.__preprocess(self.value)

        # 设置文本数据的保存路径
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = "text.json"
        # 如果路径不存在，则创建路径
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        try:
            self.text_datas[self.tag]
        except:
            self.text_datas[self.tag] = {"data": []}

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
        将文本数据写入到指定目录的指定json文件
        """
        # 将step, caption和文本数据内容写入json
        json_data = {"step": self.step, "caption": self.caption, "text": self.text_data}
        print(self.text_datas)
        self.text_datas[self.tag]["data"].append(json_data)

        try:
            with open(save_path, "w+", encoding="utf-8") as f:
                # 将文本数据写入到指定目录的指定json文件, ensure_ascii=False 保证中文不乱码
                ujson.dump(self.text_datas[self.tag], f, ensure_ascii=False)
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
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
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
