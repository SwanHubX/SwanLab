"""
@author: cunyue
@file: __init__.py
@time: 2026/3/15
@description: 视频处理模块，暂时只支持 GIF
"""

import hashlib
from io import BytesIO
from pathlib import Path
from typing import List, Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.media.video_pb2 import VideoItem, VideoValue
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg.fs import safe_write
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.video import VideoDataType

# 各格式的魔数校验表，新增格式时在此追加
# format → (magic_bytes, ...)
_FORMAT_MAGIC: dict[str, tuple[bytes, ...]] = {
    "gif": (b"GIF87a", b"GIF89a"),
}

# 路径后缀 → 格式名
_EXT_TO_FORMAT: dict[str, str] = {
    ".gif": "gif",
}


def _detect_format_by_magic(data: bytes) -> Optional[str]:
    """根据魔数推断格式，无法识别则返回 None"""
    for fmt, magics in _FORMAT_MAGIC.items():
        if any(data.startswith(m) for m in magics):
            return fmt
    return None


class Video(TransformMedia):
    def __init__(self, data_or_path: VideoDataType, caption: CaptionType = None):
        """Video class constructor

        目前支持的格式：GIF。

        Parameters
        ----------
        data_or_path: str, bytes, BytesIO, or Video
            Path to a supported video file, raw video bytes, a BytesIO containing
            video data, or another Video instance (nesting).
        caption: str, optional
            Caption for the video.
        """
        super().__init__()

        # 套娃加载
        attrs = self._unwrap(data_or_path)
        if attrs:
            self.buffer: BytesIO = attrs["buffer"]
            self.format: str = attrs["format"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        # 1. 文件路径
        if isinstance(data_or_path, str):
            ext = Path(data_or_path).suffix.lower()
            if ext not in _EXT_TO_FORMAT:
                supported = ", ".join(_EXT_TO_FORMAT)
                raise TypeError(f"Unsupported file extension '{ext}'. Supported: {supported}")
            try:
                with open(data_or_path, "rb") as f:
                    raw = f.read()
            except OSError as e:
                raise ValueError(f"Failed to open file: {data_or_path!r}") from e
            fmt = _detect_format_by_magic(raw)
            if fmt is None:
                raise TypeError(f"File '{data_or_path}' does not match any known video format magic number.")
            self.format = fmt

        # 2. bytes 或 BytesIO
        elif isinstance(data_or_path, (bytes, BytesIO)):
            raw = data_or_path if isinstance(data_or_path, bytes) else data_or_path.read()
            fmt = _detect_format_by_magic(raw)
            if fmt is None:
                supported = ", ".join(_FORMAT_MAGIC)
                raise TypeError(f"Cannot detect video format from bytes. Supported formats: {supported}")
            self.format = fmt

        # 3. 其他类型
        else:
            supported = ", ".join(_EXT_TO_FORMAT)
            raise TypeError(
                f"Unsupported type: {type(data_or_path).__name__}. "
                f"Please provide a file path ({supported}), bytes, or BytesIO."
            )

        self.buffer = BytesIO(raw)
        self.caption = caption

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_VIDEO

    @classmethod
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[VideoItem]) -> DataRecord:
        return DataRecord(
            key=key, step=step, timestamp=timestamp, type=cls.column_type(), videos=VideoValue(items=data)
        )

    def transform(self, *, step: int, path: Path) -> VideoItem:
        content = self.buffer.getvalue()
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.{self.format}"
        safe_write(path / filename, content, mode="wb")
        return VideoItem(filename=filename, sha256=sha256, size=len(content), caption=self.caption or "")
