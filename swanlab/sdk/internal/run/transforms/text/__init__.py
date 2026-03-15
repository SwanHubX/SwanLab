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

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.media.text_pb2 import TextItem, TextValue
from swanlab.sdk.internal.context import TransformMediaType
from swanlab.sdk.internal.pkg.fs import safe_write


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
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[TextItem]) -> DataRecord:
        return DataRecord(key=key, step=step, timestamp=timestamp, type=cls.column_type(), texts=TextValue(items=data))

    def transform(self, *, step: int, path: Path) -> TextItem:
        # 计算 sha256
        sha256 = hashlib.sha256(self.content.encode()).hexdigest()
        # 构建 filename
        # 历史版本直接将用户传入的content写入CH，因此这里需要增加__swanlab__后缀以区分当前版本与历史版本
        filename = f"{step:03d}-{sha256[:8]}.__swanlab__.txt"
        # 写入数据
        safe_write(path / filename, self.content)
        return TextItem(filename=filename, caption=self.caption)
