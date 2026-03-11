"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 19:17
@description: 文本处理模块
"""

import hashlib
from pathlib import Path
from typing import Optional, Union

from swanlab.proto.swanlab.data.v1.text_pb2 import TextItem, TextValue
from swanlab.sdk.internal.pkg.fs import safe_write
from swanlab.sdk.internal.run.data.transforms.abc import TransformMediaType


class Text(TransformMediaType):
    def __init__(self, content: Union["Text", str], caption: Optional[str] = None):
        attrs = self._unwrap(content)
        self.content = attrs.get("content", content)
        self.caption = caption if caption is not None else attrs.get("caption")

    @staticmethod
    def transform(
        key: str,
        step: int,
        path: Path,
        *,
        content: Union["Text", str],
        caption: Optional[str] = None,
    ) -> TextValue:
        this = Text(content=content, caption=caption)
        # 计算 sha256
        sha256 = hashlib.sha256(this.content.encode()).hexdigest()
        # 构建 filename
        # 历史版本直接将用户传入的content写入CH，因此这里需要增加__swanlab__后缀以区分当前版本与历史版本
        filename = f"{key}-{step:03d}-{sha256[:8]}.__swanlab__.txt"
        # 写入数据
        safe_write(path / filename, this.content)
        return TextValue(items=[TextItem(filename=filename, caption=this.caption)])
