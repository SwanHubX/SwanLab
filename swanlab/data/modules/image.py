import numpy as np
from PIL import Image as PILImage
from .base import BaseType
from .utils_modules import BoundingBoxes, ImageMask
from ._utils import get_file_hash_pil
from typing import Union, List, Dict, Any
from io import BytesIO
import os


def is_pytorch_tensor_typename(typename: str) -> bool:
    return typename.startswith("torch.") and ("Tensor" in typename or "Variable" in typename)


def get_full_typename(o: Any) -> Any:
    """Determine types based on type names.

    Avoids needing to to import (and therefore depend on) PyTorch, TensorFlow, etc.
    """
    instance_name = o.__class__.__module__ + "." + o.__class__.__name__
    if instance_name in ["builtins.module", "__builtin__.module"]:
        return o.__name__
    else:
        return instance_name


class Image(BaseType):
    """Image class constructor

    Parameters
    ----------
    data_or_path: (str or numpy.array or PIL.Image or torch.Tensor or matplotlib figure or List["Image"])
        Path to the image file, numpy array of image data, PIL.Image data, torch.Tensor data or matplotlib figure data.
    mode: (str)
        The PIL Mode for a image. Most commom is 'L', 'RGB', 'RGBA'.
        More information about the mode can be found at https://pillow.readthedocs.io/en/stable/handbook/concepts.html#concept-modes
    caption: (str)
        Caption for the image.
    file_type: (str)
        File type for the image. It is used to save the image in the specified format. The default is 'png'. The supported file types are ['png', 'jpg', 'jpeg', 'bmp'].
    size: (int or list or tuple)
        The size of the image can be controlled in four ways:
        1. If int type, it represents the maximum side length of the image, that is, the width and height cannot exceed this maximum side length. The image will be scaled proportionally to ensure that the maximum side length does not exceed MAX_DIMENSION.
        2. If list or tuple type with both specified values, e.g. (500, 500), then the image will be scaled to the specified width and height.
        3. If list or tuple type with only one specified value and another value as None, e.g. (500, None), it means resize the image to the specified width, and the height is scaled proportionally.
        4. If it is None, it means no scaling for the image.
    """

    def __init__(
        self,
        data_or_path: Union[str, np.ndarray, PILImage.Image, List["Image"]],
        mode: str = None,
        caption: str = None,
        file_type: str = None,
        size: Union[int, list, tuple] = None,
        # boxes: dict = None,
        # masks: dict = None,
    ):
        super().__init__(data_or_path)
        self.image_data = None
        self.mode = mode
        self.caption = self.__convert_caption(caption)
        self.format = self.__convert_file_type(file_type)
        self.size = self.__convert_size(size)

        # 提前预处理maplotlib类型
        if hasattr(self.value, "savefig"):
            # 如果输入为matplotlib图像
            self.value = self.__convert_plt_to_image(self.value)

        # TODO: 等前端支持Boxes和Masks后再开启

        # self.boxes = None
        # self.boxes_total_classes = None
        # self.masks = None
        # self.masks_total_classes = None
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
        # 判断是否要保存(mode='disabled'时不保存)
        if not self.settings.should_save:
            return
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
        return ["str", "numpy.array", "PIL.Image.Image", "torch.Tensor"]

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

    def __convert_file_type(self, file_type):
        """转换file_type，并检测file_type是否正确"""
        accepted_formats = ["png", "jpg", "jpeg", "bmp"]
        if file_type is None:
            format = "png"
        else:
            format = file_type

        if format not in accepted_formats:
            raise ValueError(f"file_type must be one of {accepted_formats}")

        return format

    def __convert_size(self, size):
        """将size转换为PIL图像的size"""
        if size is None:
            return None
        if isinstance(size, int):
            return size
        if isinstance(size, (list, tuple)):
            if len(size) == 2:
                if size[0] is None and size[1] is None:
                    return None
                elif size[0] is None:
                    return (None, int(size[1]))
                elif size[1] is None:
                    return (int(size[0]), None)
                else:
                    return (int(size[0]), int(size[1]))
            if len(size) == 1:
                if size[0] is None:
                    return None
                else:
                    return int(size[0])
        raise ValueError("size must be an int, list or tuple with 2 or 1 elements")

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
        elif is_pytorch_tensor_typename(get_full_typename(data)):
            # 如果输入为pytorch tensor
            try:
                import torchvision
            except ImportError as e:
                raise TypeError(
                    "swanlab.Image requires `torchvision` when process torch.tensor data. Install with 'pip install torchvision'."
                )

            if hasattr(data, "requires_grad") and data.requires_grad:
                data = data.detach()
            if hasattr(data, "detype") and str(data.type) == "torch.uint8":
                data = data.to(float)
            data = torchvision.utils.make_grid(data, normalize=True)
            image = PILImage.fromarray(
                data.mul(255).clamp(0, 255).byte().permute(1, 2, 0).cpu().numpy(), mode=self.mode
            )
        else:
            # 以上都不是，则报错
            raise TypeError(
                "Unsupported image type. Please provide a valid path, numpy array, PIL.Image, torch.Tensor or matplotlib figure."
            )

        self.image_data = self.__resize(image, self.size)

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
            raise TypeError(
                "Invalid numpy array for the image: the numpy array must be 2D or 3D with 3 or 4 channels."
            ) from e

    def __convert_plt_to_image(self, plt_obj):
        """ """
        try:
            buf = BytesIO()
            plt_obj.savefig(buf, format=self.format)  # 将图像保存到BytesIO对象
            buf.seek(0)  # 移动到缓冲区的开始位置
            image = PILImage.open(buf).convert(self.mode)  # 使用PIL打开图像
            buf.close()  # 关闭缓冲区
            return image
        except Exception as e:
            raise TypeError("Invalid matplotlib figure for the image") from e

    def __resize(self, image, size=None):
        """将图像调整大小"""
        # 如果size是None, 则返回原图
        if size is None:
            return image
        # 如果size是int类型，且图像的最大边长超过了size，则进行缩放
        if isinstance(size, int):
            MAX_DIMENSION = size
            if max(image.size) > MAX_DIMENSION:
                image.thumbnail((MAX_DIMENSION, MAX_DIMENSION))
        # 如果size是list或tuple类型
        elif isinstance(size, (list, tuple)):
            # 如果size是两个值的list或tuple，如(500, 500)，则进行缩放
            if None not in size:
                image = image.resize(size)
            else:
                # 如果size中有一个值为None，且图像对应的边长超过了size中的另一个值，则进行缩放
                if size[0] is not None:
                    wpercent = size[0] / float(image.size[0])
                    hsize = int(float(image.size[1]) * float(wpercent))
                    image = image.resize((size[0], hsize), PILImage.ANTIALIAS)
                elif size[1] is not None:
                    hpercent = size[1] / float(image.size[1])
                    wsize = int(float(image.size[0]) * float(hpercent))
                    image = image.resize((wsize, size[1]), PILImage.ANTIALIAS)

        return image

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
