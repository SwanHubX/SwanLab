import numpy as np
from .base import BaseType
from ..utils.file import get_text_sha256_hash
from typing import Union, List, Sequence, Tuple
import os
import ujson

NumpyHistogram = Tuple[np.ndarray, np.ndarray]


class Histogram(BaseType):
    """
    swanlab.Histogram 类型的数据，用于记录直方图数据
    """

    def __init__(
        self,
        data: Sequence = None,
        bins: int = 10,
        np_histogram: NumpyHistogram = None,
    ):
        super().__init__(data)
        self.np_histogram = np_histogram
        self.bins = bins

    def get_data(self):
        # 如果传入的是numpy直方图
        if self.np_histogram:
            if len(self.np_histogram) == 2:
                self.histogram = (
                    self.np_histogram[0].tolist() if hasattr(self.np_histogram[0], "tolist") else self.np_histogram[0]
                )
                self.bins = (
                    self.np_histogram[1].tolist() if hasattr(self.np_histogram[1], "tolist") else self.np_histogram[1]
                )
            else:
                raise TypeError(
                    "Expected np_histogram to be a tuple of (values, bin_edges) or sequence to be specified"
                )

        # 如果传入的是自定义的数据
        else:
            self.histogram, self.bins = np.histogram(
                self.value,
                bins=self.bins,
            )
            self.histogram = self.histogram.tolist()
            self.bins = self.bins.tolist()

        if len(self.histogram) + 1 != len(self.bins):
            raise TypeError("len(bins) must be len(histogram) + 1")

        self.hist_data = {"hist": self.histogram, "bins": self.bins}

        # 获取图像的hash值
        hash_name = get_text_sha256_hash(str(self.hist_data))[:16]

        # 设置保存路径, 保存文件名
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"histogram-step{self.step}-{hash_name}.json"

        # 如果不存在目录则创建
        if os.path.exists(save_dir) is False:
            os.makedirs(save_dir)
        save_path = os.path.join(save_dir, save_name)
        # 保存图像到指定目录
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["Sequence"]

    def __save(self, save_path):
        """
        将Histogram数据写入到指定目录的json文件
        """
        try:
            with open(save_path, "w", encoding="utf-8") as f:
                # 将文本数据写入到指定目录的txt文件
                writer = ujson.dump(self.hist_data, f)
        except Exception as e:
            raise TypeError(f"Could not save the table data to the path: {save_path}") from e

    def get_hist(self) -> List:
        return self.histogram

    def get_bins(self) -> List:
        return self.bins

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        return None

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Histogram"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.histogram
