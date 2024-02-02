import numpy as np
from PIL import Image as PILImage
from .base import BaseType
from typing import Union, List
import os


class Image(BaseType):
    """Image class constructor

    Parameters
    ----------
    data_or_path: (str or numpy.array or PIL.Image)
        Path to the image file, numpy array of image data or PIL.Image data.
    caption: (str)
        Caption for the image.
    """

    def __init__(self, data_or_path: Union[str, np.ndarray, PILImage.Image, List["Image"]], caption: str = None):
        super().__init__(data_or_path)
        self.image_data = None
        self.caption = caption
        if self.caption is not None:
            self.caption = self.__convert_caption(caption)

    def get_data(self):
        if isinstance(self.value, list):
            return [i.get_data() for i in self.value]
        self.__preprocess(self.value)
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"image-step{self.step}.png"
        if os.path.exists(save_dir) is False:
            os.makedirs(save_dir)
        save_path = os.path.join(save_dir, f"image-step{self.step}.png")

        # 保存图像到指定目录
        self.__save(save_path)
        return save_name

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "numpy.array", "PIL.Image.Image"]

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

    def __preprocess(self, data):
        """将不同类型的输入转换为PIL图像"""
        if isinstance(data, str):
            # 如果输入为字符串
            image = self.__load_image_from_path(data)
        elif isinstance(data, np.ndarray):
            # 如果输入为numpy array
            image = self.__convert_numpy_array_to_image(data)
        elif isinstance(data, PILImage.Image):
            # 如果输入为PIL.Image
            image = data
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported image type. Please provide a valid path, numpy array, or PIL.Image.")

        # 对数据做通道转换
        image = self.__convert_channels(image)

        # 缩放大小
        image = self.__resize(image)

        self.image_data = image

    def __load_image_from_path(self, path):
        """判断字符串是否为正确的图像路径，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.open(path)
        except Exception as e:
            raise TypeError(f"Invalid image path: {path}") from e

    def __convert_numpy_array_to_image(self, array):
        """判断np array对象是否能转换为PIL.Image，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.fromarray(array)
        except Exception as e:
            raise TypeError("Invalid numpy array for the image") from e

    def __convert_channels(self, image):
        """将1通道和4通道图像转换为3通道的RGB图像。"""
        if image.mode == "L":  # 1-channel
            return image.convert("RGB")
        elif image.mode == "RGBA" or image.mode == "CMYK":  # 4-channel
            return image.convert("RGB")
        return image

    def __resize(self, image, MAX_DIMENSION=1280):
        """将图像调整大小, 保证最大边长不超过MAX_DIMENSION"""
        if max(image.size) > MAX_DIMENSION:
            image.thumbnail((MAX_DIMENSION, MAX_DIMENSION))
        return image

    def __save(self, save_path):
        """将图像保存到指定路径"""
        pil_image = self.image_data
        if not isinstance(pil_image, PILImage.Image):
            raise TypeError("Invalid image data for the image")
        try:
            pil_image.save(save_path, format="png")
        except Exception as e:
            raise TypeError(f"Could not save the image to the path: {save_path}") from e

    def get_more(self, *args, **kwargs) -> dict:
        """返回每个step的data的更多信息的数据"""
        if isinstance(self.value, list):
            return [i.get_more() for i in self.value]
        return (
            {
                "caption": self.caption,
            }
            if self.caption is not None
            else None
        )

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Image"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.image
