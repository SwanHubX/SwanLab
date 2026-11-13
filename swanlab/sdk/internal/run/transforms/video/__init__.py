"""视频输入、同步编码与媒体记录转换。"""

from __future__ import annotations

import importlib.util
import os
from io import BytesIO
from pathlib import Path
from typing import Optional, cast

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.video import VideoErrorInfo, VideoFormat, VideoInput

from .encoding import encode, require_modules
from .source import (
    CHUNK_SIZE,
    MAX_SIZE,
    MIME_TYPES,
    VideoProcessingError,
    VideoSource,
    check_size,
    scratch_path,
    snapshot,
)


class Video(TransformMedia):
    """Log GIF/MP4/WebM/OGG files or encode TCHW/BTCHW arrays.

    ``caption`` remains the second positional argument. Arrays require an
    explicit ``format``, uint8 pixels, and use 4 fps by default. The encoded
    canvas, including batch tiling and padding, is limited to 4096 pixels per
    side. Tensor callers must detach and prepare CPU uint8 arrays themselves.
    File inputs keep their bytes and audio unchanged. Processing errors warn
    and produce an invalid object. Logging skips that video while continuing
    with other metrics.
    """

    def __init__(
        self,
        data_or_path: VideoInput,
        caption: CaptionType = None,
        fps: Optional[int] = None,
        format: Optional[VideoFormat] = None,
    ):
        super().__init__()
        self._source: Optional[VideoSource] = None
        self._error: Optional[VideoErrorInfo] = None
        self.caption: Optional[str] = None
        try:
            if caption is not None and not isinstance(caption, str):
                raise VideoProcessingError("argument", "Video caption must be a string or None")
            self.caption = caption
            self._prepare(data_or_path, fps, format)
        except Exception as error:
            self._invalidate(error)

    def _invalidate(self, error: Exception) -> None:
        code = error.code if isinstance(error, VideoProcessingError) else "encoding"
        if isinstance(error, OSError):
            code = "file_read"
        if isinstance(error, ImportError):
            code = "dependency_missing"
        self._error = VideoErrorInfo(code, str(error)[:1000])
        console.warning(f"Invalid video ({code}): {self._error.message}")

    def _prepare(self, data: object, fps: Optional[int], fmt: Optional[VideoFormat]) -> None:
        if isinstance(data, Video):
            self._source = data._source
            self._error = data.error
            if self.caption is None:
                self.caption = data.caption
            self._file_options(fps, fmt, data.format)
            return
        if isinstance(data, (str, os.PathLike)):
            filename = os.fspath(data)
            if not isinstance(filename, str):
                raise VideoProcessingError("input_type", "Video paths must contain text")
            path = Path(filename)
            actual_format = path.suffix.lower().lstrip(".")
            if actual_format not in MIME_TYPES:
                raise VideoProcessingError("format_mismatch", "Video file extension must be gif, mp4, webm or ogg")
            check_size(path.stat().st_size)
            self._file_options(fps, fmt, actual_format)
            with path.open("rb") as stream:
                self._source = snapshot(stream, actual_format, scratch_path(actual_format))
            return
        if isinstance(data, (bytes, BytesIO)):
            if isinstance(data, bytes):
                check_size(len(data))
            stream = BytesIO(data) if isinstance(data, bytes) else data
            # 不移动调用方游标来探测头部；getbuffer 的切片不复制整个输入。
            view = stream.getbuffer()
            try:
                is_gif = bytes(view[stream.tell() : stream.tell() + 6]) in (b"GIF87a", b"GIF89a")
            finally:
                view.release()
            actual_format = fmt or ("gif" if is_gif else None)
            if actual_format is None:
                consumed = 0
                while chunk := stream.read(min(CHUNK_SIZE, MAX_SIZE + 1 - consumed)):
                    consumed += len(chunk)
                    check_size(consumed)
                raise VideoProcessingError("format_required", "Non-GIF video bytes require an explicit format")
            if actual_format not in MIME_TYPES:
                raise VideoProcessingError("argument", "Video format must be gif, mp4, webm or ogg")
            self._file_options(fps, None, actual_format)
            self._source = snapshot(stream, actual_format, scratch_path(actual_format))
            return
        if fmt not in MIME_TYPES:
            raise VideoProcessingError(
                "format_required", "Video arrays require an explicit gif, mp4, webm or ogg format"
            )
        converter = getattr(data, "numpy", None)
        if callable(converter):
            try:
                data = converter()
            except Exception as error:
                raise VideoProcessingError(
                    "numpy_conversion",
                    f"Video numpy() conversion failed: {type(error).__name__}; "
                    "provide a CPU uint8 NumPy array, detaching tensors first",
                ) from error
        elif importlib.util.find_spec("numpy") is None:
            raise VideoProcessingError("input_type", f"Unsupported video input: {type(data).__name__}")
        require_modules("numpy")
        import numpy as np

        if not isinstance(data, np.ndarray):
            raise VideoProcessingError("input_type", "Video requires a path, bytes, ndarray or numpy() object")
        self._source = encode(data, fmt, fps, scratch_path(fmt))

    @staticmethod
    def _file_options(fps: Optional[int], fmt: Optional[str], actual: Optional[str]) -> None:
        if fps is not None:
            console.warning("Video fps is ignored for existing video content")
        if fmt is not None and fmt != actual:
            console.warning("Video format is determined by the existing file or nested video; format argument ignored")

    @property
    def is_valid(self) -> bool:
        return self._source is not None and self._error is None

    @property
    def error(self) -> Optional[VideoErrorInfo]:
        return self._error

    @property
    def format(self) -> Optional[VideoFormat]:
        return cast(Optional[VideoFormat], self._source.format if self._source is not None else None)

    @property
    def buffer(self) -> BytesIO:
        if not self.is_valid:
            return BytesIO()
        assert self._source is not None
        try:
            return self._source.buffer
        except Exception as error:
            self._invalidate(error)
            return BytesIO()

    @buffer.setter
    def buffer(self, value: BytesIO) -> None:
        try:
            if not isinstance(value, BytesIO) or self._source is None:
                raise VideoProcessingError("argument", "Video buffer must be a BytesIO for an existing video source")
            source = VideoSource(self._source.format, buffer=value)
            source.refresh()
            # 整体替换仅作用于当前包装对象，已有包装对象继续持有原数据源。
            self._source = source
            self._error = None
        except Exception as error:
            self._invalidate(error)

    @property
    def size(self) -> Optional[int]:
        return self._source.size if self.is_valid and self._source is not None else None

    @property
    def sha256(self) -> Optional[str]:
        return self._source.sha256 if self.is_valid and self._source is not None else None

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_VIDEO

    def transform(self, *, step: int, path: Optional[Path]) -> MediaItem:
        """Build a media item; the video record builder isolates conversion failures."""
        if not self.is_valid:
            error = self._error
            raise VideoProcessingError(
                error.code if error is not None else "invalid_video",
                error.message if error is not None else "Video has no valid source",
            )
        assert self._source is not None
        try:
            self._source.refresh()
            item = MediaItem(
                filename=f"{step:03d}-{self._source.sha256}.{self.format}",
                sha256=self._source.sha256,
                size=self._source.size,
                caption=self.caption or "",
            )
            if path is None:
                view = self._source.buffer.getbuffer()
                try:
                    item.payload = bytes(view)
                finally:
                    view.release()
            else:
                self._source.store(path / item.filename)
            return item
        except Exception as error:
            self._invalidate(error)
            raise
