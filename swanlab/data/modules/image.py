import numpy as np
from PIL import Image as PILImage
from .base import BaseType
from .utils_modules import BoundingBoxes, ImageMask
from ..utils.file import get_file_hash_pil
from typing import Union, List, Dict
from io import BytesIO
import os


class Image(BaseType):
    """Image class constructor

    Parameters
    ----------
    data_or_path: (str or numpy.array or PIL.Image)
        Path to the image file, numpy array of image data or PIL.Image data.
    mode: (str)
        The PIL Mode for a image. Most commom is 'L', 'RGB', 'RGBA'.
        More information about the mode can be found at https://pillow.readthedocs.io/en/stable/handbook/concepts.html#concept-modes
    caption: (str)
        Caption for the image.
    """

    def __init__(
        self,
        data_or_path: Union[str, np.ndarray, PILImage.Image, List["Image"]],
        mode: str = "RGB",
        caption: str = None,
        # boxes: dict = None,
        # masks: dict = None,
    ):
        super().__init__(data_or_path)
        self.image_data = None
        self.mode = mode
        self.caption = self.__convert_caption(caption)

        # self.boxes = None
        # self.boxes_total_classes = None
        # self.masks = None
        # self.masks_total_classes = None

        # TODO: 等前端支持Boxes和Masks后再开启
        # if boxes:
        #     self.boxes, self.boxes_total_classes = self.__convert_boxes(boxes)

        # if masks:
        #     self.masks, self.masks_total_classes = self.__convert_masks(masks)

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
        save_name = (
            f"{self.caption}-step{self.step}-{hash_name}.png"
            if self.caption is not None
            else f"image-step{self.step}-{hash_name}.png"
        )
        # 如果不存在目录则创建
        if os.path.exists(save_dir) is False:
            os.makedirs(save_dir)
        save_path = os.path.join(save_dir, save_name)
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
        # 如果类型是None，则转换为默认字符串
        elif caption is None:
            caption = None
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption

    def __convert_boxes(self, boxes):
        """将boxes转换为Dict[str, BoundingBoxes]类型对象, 并返回该对象和总标签"""
        # 如果boxes的类型不是字典，则报错
        if not isinstance(boxes, dict):
            raise TypeError("swanlab.Image 'boxes' argument must be a dictionary")

        boxes_final: Dict[str, BoundingBoxes] = {}
        total_classes = {}

        # 对于boxes中的每一个key
        for key in boxes:
            box_item = boxes[key]
            if isinstance(box_item, BoundingBoxes):
                boxes_final[key] = box_item
            elif isinstance(box_item, dict):
                boxes_final[key] = BoundingBoxes(box_item, key)
            total_classes.update(boxes_final[key]._class_labels)

        return boxes_final, total_classes

    def __convert_masks(self, masks):
        """将masks转换为Dict[str, ImageMask]类型对象, 并返回该对象和总标签"""
        # 如果masks的类型不是字典，则报错
        if not isinstance(masks, dict):
            raise TypeError("swanlab.Image 'masks' argument must be a dictionary")

        masks_final: Dict[str, ImageMask] = {}
        total_classes = {}

        # 对于masks中的每一个key
        for key in masks:
            mask_item = masks[key]
            if isinstance(mask_item, ImageMask):
                masks_final[key] = mask_item
            elif isinstance(mask_item, dict):
                masks_final[key] = ImageMask(mask_item, key)
            total_classes.update(masks_final[key]._class_labels)

        return masks_final, total_classes

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
            image = data.convert(self.mode)
        elif hasattr(self.value, "savefig"):
            # 如果输入为matplotlib图像
            image = self.__convert_plt_to_image(data)
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported image type. Please provide a valid path, numpy array, or PIL.Image.")
        # 缩放大小
        image = self.__resize(image)

        self.image_data = image

    def __load_image_from_path(self, path):
        """判断字符串是否为正确的图像路径，如果是则返回np.ndarray类型对象，如果不是则报错"""
        try:
            return PILImage.open(path).convert(self.mode)
        except Exception as e:
            raise TypeError(f"Invalid image path: {path}") from e

    def __convert_numpy_array_to_image(self, array):
        """ """
        try:
            if array.ndim == 2 or (array.ndim == 3 and array.shape[2] in [3, 4]):
                array = np.clip(array, 0, 255).astype(np.uint8)
                return PILImage.fromarray(array, mode=self.mode)
            else:
                raise TypeError("Invalid numpy array: the numpy array must be 2D or 3D with 3 or 4 channels.")
        except Exception as e:
            raise TypeError("Invalid numpy array for the image") from e

    def __convert_plt_to_image(self, plt_obj):
        """ """
        try:
            buf = BytesIO()
            plt_obj.savefig(buf, format="png")  # 将图像保存到BytesIO对象
            buf.seek(0)  # 移动到缓冲区的开始位置
            image = PILImage.open(buf).convert(self.mode)  # 使用PIL打开图像
            buf.close()  # 关闭缓冲区
            return image
        except Exception as e:
            raise TypeError("Invalid matplotlib figure for the image") from e

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
        """返回config数据"""
        # 如果传入的是Image类列表
        if isinstance(self.value, list):
            return self.get_more_list()
        else:
            get_more_dict = {}

            if self.caption is not None:
                get_more_dict["caption"] = self.caption

            # if self.boxes is not None:
            #     get_more_dict["boxes"] = self.boxes

            # if self.boxes_total_classes is not None:
            #     get_more_dict["boxes_total_classes"] = self.boxes_total_classes

            # if self.masks is not None:
            #     get_more_dict["masks"] = self.masks

            # if self.masks_total_classes is not None:
            #     get_more_dict["masks_total_classes"] = self.masks_total_classes

            return get_more_dict if get_more_dict else None

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Image"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.image
