import numpy as np
from PIL import Image as PILImage
from .base import BaseType
import os


class Image(BaseType):
    def __init__(self, data_or_path):
        super().__init__(data_or_path)
        self.image_data = None

    def get_data(self):
        self.preprocess(self.value)
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"image-step{self.step}.png"
        if os.path.exists(self.settings.static_dir) is False:
            os.makedirs(self.settings.static_dir)
        save_path = os.path.join(save_dir, f"image-step{self.step}.png")

        # 保存图像到指定目录
        self.save(save_path)
        return save_name

    def preprocess(self, data):
        """将不同类型的输入转换为PIL图像"""
        if isinstance(data, str):
            # 如果输入为字符串
            image = self.load_image_from_path(data)
        elif isinstance(data, np.ndarray):
            # 如果输入为numpy array
            image = self.convert_numpy_array_to_image(data)
        elif isinstance(data, PILImage.Image):
            # 如果输入为PIL.Image
            image = data
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported image type. Please provide a valid path, numpy array, or PIL.Image.")

        # 对数据做通道转换
        image = self.convert_channels(image)

        # 缩放大小
        image = self.resize(image)

        self.image_data = image

    def load_image_from_path(self, path):
        """判断字符串是否为正确的图像路径，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.open(path)
        except Exception as e:
            raise ValueError(f"Invalid image path: {path}") from e

    def convert_numpy_array_to_image(self, array):
        """判断np array对象是否能转换为PIL.Image，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.fromarray(array)
        except Exception as e:
            raise ValueError("Invalid numpy array for the image") from e

    def convert_channels(self, image):
        """将1通道和4通道图像转换为3通道的RGB图像。"""
        if image.mode == "L":  # 1-channel
            return image.convert("RGB")
        elif image.mode == "RGBA" or image.mode == "CMYK":  # 4-channel
            return image.convert("RGB")
        return image

    def resize(self, image, MAX_DIMENSION=1280):
        """将图像调整大小, 保证最大边长不超过MAX_DIMENSION"""
        if max(image.size) > MAX_DIMENSION:
            image.thumbnail((MAX_DIMENSION, MAX_DIMENSION))
        return image

    def save(self, save_path):
        """将图像保存到指定路径"""
        try:
            self.image.save(save_path)
        except Exception as e:
            raise ValueError(f"Could not save the image to the path: {save_path}") from e

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Image"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.image
