import numpy as np
from PIL import Image as PILImage
from .base import BaseType
from ..utils.file import get_file_hash_pil
from typing import Union, List, Sequence, Tuple
import os

NumpyHistogram = Tuple[np.ndarray, np.ndarray]


class Histogram(BaseType):
    """
    swanlab.Histogram 类型的数据，用于记录直方图数据

    """

    def __init__(
        self,
        data: Union[Sequence, List["Histogram"]],
        bins: int = 10,
        np_histogram: NumpyHistogram = None,
    ):

        super().__init__(data)
        # 如果传入的是numpy直方图
        if np_histogram:
            if len(np_histogram) == 2:
                self.histogram = np_histogram[0].tolist() if hasattr(np_histogram[0], "tolist") else np_histogram[0]
                self.bins = np_histogram[1].tolist() if hasattr(np_histogram[1], "tolist") else np_histogram[1]
            else:
                raise ValueError(
                    "Expected np_histogram to be a tuple of (values, bin_edges) or sequence to be specified"
                )
        # 如果传入的是自定义的数据
        else:
            self.histogram, self.bins = np.histogram(
                data,
                bins=bins,
            )
            self.histogram = self.histogram.tolist()
            self.bins = self.bins.tolist()

        if len(self.histogram) + 1 != len(self.bins):
            raise ValueError("len(bins) must be len(histogram) + 1")

    def get_data(self):
        # 如果传入的是Image类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        # 图像预处理
        self.__preprocess(self.value)
        # 获取图像的hash值
        hash_name = get_file_hash_pil(self.image_data)[:16]
        # 设置保存路径, 保存文件名
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"image-step{self.step}-{hash_name}.{self.format}"
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
        """将图像保存到指定路径"""
        pil_image = self.image_data
        if not isinstance(pil_image, PILImage.Image):
            raise TypeError("Invalid image data for the image")
        try:
            if self.format == "jpg":
                pil_image.save(save_path, format="JPEG")
            else:
                pil_image.save(save_path, format=self.format)

        except Exception as e:
            raise TypeError(f"Could not save the image to the path: {save_path}") from e

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        return None

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Histogram"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.histogram
