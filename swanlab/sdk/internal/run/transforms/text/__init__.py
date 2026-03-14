"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 19:17
@description: 文本处理模块
"""

import hashlib
from pathlib import Path
from typing import List, Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.data.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.data.v1.metric_pb2 import MetricRecord
from swanlab.proto.swanlab.data.v1.text_pb2 import TextItem, TextValue
from swanlab.sdk.internal.context import TransformMediaType
from swanlab.sdk.internal.pkg.fs import safe_write
from swanlab.sdk.typings.run.data import MediaTransferType


class Text(TransformMediaType):
    def __init__(self, content: Union["Text", str], caption: Optional[str] = None):
        super().__init__()
        attrs = self._unwrap(content)
        self.content = attrs.get("content", content)
        self.caption = caption if caption is not None else attrs.get("caption")

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_TEXT

    @classmethod
    def type(cls) -> MediaTransferType:
        return "text"

    @classmethod
    def build_metric_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[TextItem]) -> MetricRecord:
        return MetricRecord(key=key, step=step, timestamp=timestamp, texts=TextValue(items=data, path="text"))

    def transform(self, key: str, step: int, path: Path) -> TextItem:
        # 计算 sha256
        sha256 = hashlib.sha256(self.content.encode()).hexdigest()
        # 构建 filename
        # 历史版本直接将用户传入的content写入CH，因此这里需要增加__swanlab__后缀以区分当前版本与历史版本
        filename = f"{key}-{step:03d}-{sha256[:8]}.__swanlab__.txt"
        # 写入数据
        safe_write(path / filename, self.content)
        return TextItem(filename=filename, caption=self.caption)
