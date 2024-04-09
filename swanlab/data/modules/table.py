# -*- coding: utf-8 -*-
from .base import BaseType
import os
from typing import List
import ujson
from ..utils.file import get_text_sha256_hash


class Table(BaseType):
    """Table class constructor

    Parameters
    ----------
    columns: List[str]
        Names of the columns in the Table.
    data: List[List[str]]
        2D row-oriented array of values.
    dataframe: pd.DataFrame
        DataFrame object used to create the Table. When set, `data` and `columns` arguments are ignored.
    """

    def __init__(self, columns: List[str] = None, data: List[List[str]] = None, dataframe=None):

        # 如果dataframe存在，那么将无视columns和data, 从dataframe中提取数据
        if dataframe is not None:
            columns, data = self.__init_from_pd_dataframe(dataframe)

        super().__init__(data)

        # 对colums做类型检查
        if not isinstance(columns, list):
            raise TypeError("columns must be a list.")
        if not all(isinstance(column, str) for column in columns):
            raise TypeError("columns must be a list of strings.")
        if len(columns) == 0:
            raise TypeError("columns must not be empty.")
        self.columns = columns

        self.colums_length = len(self.columns)
        self.tb_data = None

    def get_data(self):
        # 预处理文本数据
        self.__preprocess(self.value)

        # 设置文本数据的保存路径
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        # 将tb_data转换为字符串
        tb_data_str = "".join([str(item) for sublist in self.tb_data for item in sublist])

        hash_name = get_text_sha256_hash(tb_data_str)[:16]
        save_name = f"table-{self.tag}-{hash_name}.json"

        # 如果路径不存在，则创建路径
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存文本数据写入到指定目录的指定json文件
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["list"]

    def __init_from_pd_dataframe(self, dataframe):
        try:
            import pandas as pd
        except ImportError as e:
            raise TypeError("swanlab.Table requires pandas when process dataframe. Install with 'pip install pandas'.")

        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError("data_frame must be a pandas DataFrame.")
        # 提取列名作为表头
        columns = dataframe.columns.tolist()

        # 将DataFrame的行转换为列表的列表（二维列表）
        data = dataframe.values.tolist()

        return columns, data

    def __preprocess(self, data: List[List[str]]):
        """
        根据不同的输入类型进行不同处理
        """
        # 必须是列表
        if not isinstance(data, list):
            raise TypeError("data must be a list.")
        # 必须不为空
        if len(data) == 0:
            raise TypeError("data must not be empty.")
        # 必须是二维数组
        if not all(isinstance(row, list) for row in data):
            raise TypeError("data must be a 2D row-oriented array.")

        data_list_convert_string = []
        # 遍历二维数组，将元素转换为字符串
        for row in data:
            # 检查每一行的长度是否等于columns的长度
            if len(row) != self.colums_length:
                raise TypeError("The length of data's row must be equal to the length of columns.")
            try:
                # 将每个元素转换为字符串
                row = [str(item) for item in row]
                data_list_convert_string.append(row)
            except Exception as e:
                raise TypeError(
                    "Elements within 'data' in swanlab.Table contain values that cannot be converted to strings."
                ) from e

        self.tb_data = data_list_convert_string

    def __save(self, save_path):
        """
        将Table数据写入到指定目录的json文件
        """
        data = {"columns": self.columns, "data": self.tb_data}

        try:
            with open(save_path, "w", encoding="utf-8") as f:
                # 将文本数据写入到指定目录的txt文件
                writer = ujson.dump(data, f)
        except Exception as e:
            raise TypeError(f"Could not save the table data to the path: {save_path}") from e

    def get_columns(self) -> list:
        """返回列名"""
        return self.columns

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        return None

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Table"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.table
