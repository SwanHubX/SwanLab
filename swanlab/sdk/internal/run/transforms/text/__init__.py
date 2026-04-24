"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 19:17
@description: 文本处理模块
"""

import hashlib
from pathlib import Path

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.text import TextDataType


class Text(TransformMedia):
    def __init__(self, content: TextDataType, caption: CaptionType = None):
        super().__init__()
        attrs = self._unwrap(content)
        self.content = attrs.get("content", content)
        self.caption = caption if caption is not None else attrs.get("caption")

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_TEXT

    def transform(self, *, step: int, path: Path) -> MediaItem:
        content_encode = self.content.encode()
        # 计算 sha256
        sha256 = hashlib.sha256(content_encode).hexdigest()
        # 构建 filename
        # 历史版本直接将用户传入的content写入CH，这交给前端去适配
        filename = f"{step:03d}-{sha256[:8]}.txt"
        # 写入数据
        fs.safe_write(path / filename, self.content)
        return MediaItem(filename=filename, sha256=sha256, size=len(content_encode), caption=self.caption)
