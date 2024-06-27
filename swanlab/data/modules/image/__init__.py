from swankit.core import MediaBuffer, DataSuite as D, MediaType
from typing import Union, Any, TYPE_CHECKING
from io import BytesIO

if TYPE_CHECKING:
    # noinspection PyPackageRequirements
    import matplotlib as _matplotlib  # type: ignore

    # noinspection PyPackageRequirements
    import numpy as _numpy  # type: ignore

    # noinspection PyPackageRequirements
    import torch as _torch  # type: ignore

    # noinspection PyPackageRequirements
    from PIL.Image import Image as _PILImage  # type: ignore

    TorchTensorType = Union["_torch.Tensor", "_torch.Variable"]
    ImageDataType = Union["_matplotlib.artist.Artist", "_PILImage", "TorchTensorType", "_numpy.ndarray"]
    ImageDataOrPathType = Union[str, "_PILImage", ImageDataType]


def is_pytorch_tensor_typename(typename: str) -> bool:
    return typename.startswith("torch.") and ("Tensor" in typename or "Variable" in typename)


def get_full_typename(o: Any) -> Any:
    """Determine types based on type names.

    Avoids needing to import (and therefore depend on) PyTorch, TensorFlow, etc.
    """
    instance_name = o.__class__.__module__ + "." + o.__class__.__name__
    if instance_name in ["builtins.module", "__builtin__.module"]:
        return o.__name__
    else:
        return instance_name


def convert_size(size=None):
    """将size转换为PIL图像的size"""
    if size is None:
        return None

    elif isinstance(size, int):
        return size

    elif isinstance(size, (list, tuple)) and len(size) in [1, 2]:
        size = tuple(size)
        width = int(size[0]) if size[0] is not None else None
        height = int(size[1]) if len(size) > 1 and size[1] is not None else None

        return (width, height) if len(size) == 2 else width

    raise ValueError("swanlab.Image - param `size` must be a list (with 2 or 1 elements) or an int")


class Image(MediaType):
    ACCEPT_FORMAT = ["png", "jpg", "jpeg", "bmp"]

    def __init__(
        self,
        data_or_path: "ImageDataOrPathType",
        mode: str = None,
        caption: str = None,
        file_type: str = None,
        size: Union[int, list, tuple] = None,
    ):
        """Image class constructor

        Parameters
        ----------
        data_or_path: (str or numpy.ndarray or PIL.Image or torch.Tensor or matplotlib figure)
            Path to the image file, numpy array of image data, PIL.Image data, torch.Tensor data or
            matplotlib figure data.
        mode: (str)
            The PIL Mode for an image. Most common is 'L', 'RGB', 'RGBA'.
            More information about the mode can be found at
            https://pillow.readthedocs.io/en/stable/handbook/concepts.html#concept-modes
        caption: (str)
            Caption for the image.
        file_type: (str)
            File type for the image. It is used to save the image in the specified format. The default is 'png'.
            The supported file types are ['png', 'jpg', 'jpeg', 'bmp'].
        size: (int or list or tuple)
            The size of the image can be controlled in four ways:
            If int type, it represents the maximum side length of the image, that is,
            the width and height cannot exceed this maximum side length. The image will be scaled proportionally to
            ensure that the maximum side length does not exceed MAX_DIMENSION;

            If list or tuple type with both specified values, e.g. (500, 500),
            then the image will be scaled to the specified width and height;

            If list or tuple type with only one specified value and another value as None, e.g. (500, None),
            it means resize the image to the specified width, and the height is scaled proportionally.

            If it is None, it means no scaling for the image.
        """

        try:
            # noinspection PyPackageRequirements
            from PIL import Image as PILImage

            # noinspection PyPackageRequirements
            import numpy as np
        except ImportError:
            raise ImportError(
                "pillow and numpy are required for Image class, you can install them by `pip install pillow numpy`"
            )

        super().__init__()
        self.format = self.__convert_file_type(file_type)
        self.size = convert_size(size)

        if isinstance(data_or_path, str):
            # 如果输入为路径字符串
            try:
                image_data = PILImage.open(data_or_path).convert(mode)
            except Exception as e:
                raise ValueError(f"Invalid image path: {data_or_path}") from e
        elif isinstance(data_or_path, PILImage.Image):
            # 如果输入为PIL.Image
            image_data = data_or_path.convert(mode)
        elif is_pytorch_tensor_typename(get_full_typename(data_or_path)):
            # 如果输入为pytorch tensor
            try:
                import torchvision  # noqa
            except ImportError:
                raise TypeError(
                    "swanlab.Image requires `torchvision` when process torch.tensor data. "
                    "Install with 'pip install torchvision'."
                )
            if hasattr(data_or_path, "requires_grad") and data_or_path.requires_grad:
                data_or_path = data_or_path.detach()  # noqa
            if hasattr(data_or_path, "detype") and str(data_or_path.type) == "torch.uint8":
                data_or_path = data_or_path.to(float)  # type: ignore
            data_or_path = torchvision.utils.make_grid(data_or_path, normalize=True)
            image_data = PILImage.fromarray(
                data_or_path.mul(255).clamp(0, 255).byte().permute(1, 2, 0).cpu().numpy(), mode=mode
            )
        elif hasattr(data_or_path, "savefig"):
            # 如果输入为matplotlib图像
            try:
                # noinspection PyPackageRequirements
                import matplotlib
            except ImportError:
                raise ImportError(
                    "swanlab.Image requires `matplotlib` when process matplotlib.artist.Artist data. "
                    "you can install them by `pip install matplotlib`"
                )

            try:
                buf = BytesIO()
                data_or_path.savefig(buf, format=self.format)  # 将图像保存到BytesIO对象
                buf.seek(0)  # 移动到缓冲区的开始位置
                image_data = PILImage.open(buf).convert(mode)  # 使用PIL打开图像
                buf.close()  # 关闭缓冲区
            except Exception as e:
                raise TypeError("Invalid matplotlib figure for the image") from e

        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy array
            try:
                if data_or_path.ndim == 2 or (data_or_path.ndim == 3 and data_or_path.shape[2] in [3, 4]):
                    image_data = PILImage.fromarray(np.clip(data_or_path, 0, 255).astype(np.uint8), mode=mode)
                else:
                    raise TypeError("Invalid numpy array: the numpy array must be 2D or 3D with 3 or 4 channels.")
            except Exception as e:
                raise TypeError("Invalid numpy array for the image") from e

        else:
            # 以上都不是，则报错
            raise TypeError(
                "Unsupported image type. Please provide a valid path, numpy array, PIL.Image, torch."
                "Tensor or matplotlib figure."
            )

        self.image_data = self.__resize(image_data, self.size)
        """
        转换为矩阵后的数据
        """
        self.buffer = MediaBuffer()
        self.image_data.save(self.buffer, format=self.format if self.format != "jpg" else "jpeg")
        self.caption = D.check_caption(caption)

    def __convert_file_type(self, file_type: str = None):
        """转换file_type，并检测file_type是否正确"""
        if file_type is None:
            file_type = "png"

        if file_type not in self.ACCEPT_FORMAT:
            raise ValueError(f"file_type must be one of {self.ACCEPT_FORMAT}")

        return file_type

    def __resize(self, image, size=None):
        """将图像调整大小"""
        if size is None:
            self.size = image.size
            return image

        elif isinstance(size, int):
            if max(image.size) > size:
                image.thumbnail((size, size))
            return image

        elif isinstance(size, (list, tuple)):
            size = tuple(size)
            if all(size):
                return image.resize(size)
            if size[0] is not None:
                w_percent = size[0] / float(image.size[0])
                hsize = int(float(image.size[1]) * w_percent)
                return image.resize((size[0], hsize))
            if size[1] is not None:
                h_percent = size[1] / float(image.size[1])
                w_size = int(float(image.size[0]) * h_percent)
                return image.resize((w_size, size[1]))

        raise ValueError("swanlab.Image - param `size` must be a list (with 2 or 1 elements) or an int")

    @property
    def image_size(self):
        """
        得到图像的大小
        """
        return self.image_data.size

    # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        hash_name = D.get_hash_by_pil(self.image_data)[:16]
        save_name = f"image-step{self.step}-{hash_name}.{self.format}"
        return save_name, self.buffer

    def get_more(self):
        """返回more数据"""
        return {"caption": self.caption} if self.caption is not None else None

    def get_section(self) -> str:
        """设定分组名"""
        return "Image"

    def get_chart(self):
        return self.Chart.IMAGE
