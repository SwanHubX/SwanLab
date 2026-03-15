"""
@author: cunyue
@file: __init__.py
@time: 2026/3/15
@description: 图像处理模块
"""

import hashlib
from io import BytesIO
from pathlib import Path
from typing import List, Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.media.image_pb2 import ImageItem, ImageValue
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg.fs import safe_write

ACCEPT_FORMAT = ["png", "jpg", "jpeg", "bmp"]


def _is_torch_tensor(obj) -> bool:
    """通过类型名检测 PyTorch Tensor，避免强制导入 torch"""
    typename = obj.__class__.__module__ + "." + obj.__class__.__name__
    return typename.startswith("torch.") and ("Tensor" in typename or "Variable" in typename)


def _resize(image: "vendor.PIL.Image.Image", size) -> "vendor.PIL.Image.Image":
    """按 size 参数缩放图像"""
    if size is None:
        return image
    if isinstance(size, int):
        if max(image.size) > size:
            image.thumbnail((size, size))
        return image
    if isinstance(size, (list, tuple)):
        w, h = (tuple(size) + (None,))[:2]
        if w is not None and h is not None:
            return image.resize((int(w), int(h)))
        if w is not None:
            return image.resize((int(w), int(image.size[1] * w / image.size[0])))
        if h is not None:
            return image.resize((int(image.size[0] * h / image.size[1]), int(h)))
    raise ValueError("size must be an int, or a list/tuple with 1-2 elements")


class Image(TransformMedia):
    def __init__(
        self,
        data_or_path: Union[
            "Image",
            str,
            "vendor.PIL.Image.Image",
            "vendor.np.ndarray",
            "vendor.torch.Tensor",
            "vendor.matplotlib.figure.Figure",
        ],
        mode: Optional[str] = None,
        caption: Optional[str] = None,
        file_type: Optional[str] = None,
        size: Optional[Union[int, list, tuple]] = None,
    ):
        """Image class constructor

        Parameters
        ----------
        data_or_path: str, PIL.Image.Image, numpy.ndarray, torch.Tensor, matplotlib.figure.Figure, or Image
            Path to an image file (PNG/JPG/JPEG/BMP; GIF is not supported), a PIL Image,
            numpy array (shape: (H, W) or (H, W, 3/4)), torch.Tensor, matplotlib figure,
            or another Image instance (nesting).
        mode: str, optional
            PIL mode applied when converting to PIL.Image (e.g. 'RGB', 'L').
        caption: str, optional
            Caption for the image.
        file_type: str, optional
            Output file format. One of ['png', 'jpg', 'jpeg', 'bmp']. Defaults to 'png'.
        size: int, list, or tuple, optional
            Resize policy:
            - int: maximum side length (aspect-ratio preserved via thumbnail).
            - (w, h): exact target size.
            - (w, None) / (None, h): fix one dimension, scale the other proportionally.
            - None: no resize.
        """
        super().__init__()

        # 套娃加载
        attrs = self._unwrap(data_or_path)
        if attrs:
            self.buffer: BytesIO = attrs["buffer"]
            self.file_type: str = attrs["file_type"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        # 校验 file_type
        ft = (file_type or "png").lower()
        if ft not in ACCEPT_FORMAT:
            raise ValueError(f"Unsupported file_type '{ft}'. Accepted: {ACCEPT_FORMAT}")
        self.file_type = ft

        # ---------- 各类型输入 → PIL Image ----------
        # 考虑到懒加载限制我们一般使用鸭子类型判断，而不是 isinstance
        PILImage = vendor.PIL.Image
        # 1. 文件路径
        if isinstance(data_or_path, str):
            if data_or_path.lower().endswith(".gif"):
                raise TypeError("GIF images are not supported. Please convert to PNG or JPG first.")
            try:
                pil_img = PILImage.open(data_or_path)
            except Exception as e:
                raise ValueError(f"Failed to open image file: {data_or_path!r}") from e
            if getattr(pil_img, "format", None) == "GIF":
                raise TypeError("GIF images are not supported. Please convert to PNG or JPG first.")
            image_data = pil_img.convert(mode)
        # 2. PIL Image
        elif isinstance(data_or_path, PILImage.Image):
            if getattr(data_or_path, "format", None) == "GIF":
                raise TypeError("GIF images are not supported. Please convert to PNG or JPG first.")
            image_data = data_or_path.convert(mode)
        # 3. PyTorch Tensor
        elif _is_torch_tensor(data_or_path):
            t = data_or_path
            if hasattr(t, "requires_grad") and t.requires_grad:  # type: ignore[union-attr]
                t = t.detach()  # type: ignore[union-attr]
            t = vendor.torchvision.utils.make_grid(t, normalize=True)  # type: ignore[arg-type]
            image_data = PILImage.fromarray(t.mul(255).clamp(0, 255).byte().permute(1, 2, 0).cpu().numpy(), mode=mode)
        # 4. Matplotlib Figure
        elif hasattr(data_or_path, "savefig"):
            try:
                buf = BytesIO()
                data_or_path.savefig(buf, format=self.file_type)  # type: ignore[union-attr]
                buf.seek(0)
                image_data = PILImage.open(buf).convert(mode)
                buf.close()
            except Exception as e:
                raise TypeError("Failed to convert matplotlib figure to image") from e
        # 5. Numpy Array
        elif isinstance(data_or_path, vendor.np.ndarray):
            arr = data_or_path
            if arr.ndim == 2 or (arr.ndim == 3 and arr.shape[2] in (3, 4)):
                image_data = PILImage.fromarray(vendor.np.clip(arr, 0, 255).astype(vendor.np.uint8), mode=mode)
            else:
                raise TypeError(f"Invalid numpy array shape for Image: expected (H, W) or (H, W, 3/4), got {arr.shape}")

        else:
            raise TypeError(
                f"Unsupported image type: {type(data_or_path).__name__}. "
                "Please provide a file path, PIL.Image, numpy.ndarray, torch.Tensor, or matplotlib figure."
            )

        image_data = _resize(image_data, size)
        self.buffer = BytesIO()
        save_fmt = "jpeg" if self.file_type == "jpg" else self.file_type
        image_data.save(self.buffer, format=save_fmt)
        self.caption = caption

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_IMAGE

    @classmethod
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[ImageItem]) -> DataRecord:
        return DataRecord(
            key=key, step=step, timestamp=timestamp, type=cls.column_type(), images=ImageValue(items=data)
        )

    def transform(self, *, step: int, path: Path) -> ImageItem:
        content = self.buffer.getvalue()
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.{self.file_type}"
        safe_write(path / filename, content, mode="wb")
        return ImageItem(filename=filename, sha256=sha256, size=len(content), caption=self.caption or "")
