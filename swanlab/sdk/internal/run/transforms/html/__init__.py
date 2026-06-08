"""
@author: caddiesnew
@file: __init__.py
@time: 2026/6/8
@description: HTML 数据处理模块
"""

import hashlib
from pathlib import Path
from typing import Optional

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.html import HtmlDataType


class Html(TransformMedia):
    def __init__(self, data: HtmlDataType, caption: CaptionType = None):
        """
        HTML 数据类，支持多种输入格式。

        Parameters
        ----------
        data: str, pathlib.Path, TextIO, or Html
            - str: 可以是 HTML 文件路径（.html 后缀且文件存在）或原始 HTML 字符串
            - pathlib.Path: HTML 文件路径
            - TextIO: 文件类对象，如 open("file.html") 返回的对象
            - Html: 套娃加载，从另一个 Html 实例复制属性
        caption: str, optional
            可选的标题说明
        """
        super().__init__()

        # 套娃加载
        attrs = self._unwrap(data)
        if attrs:
            self.content: str = attrs["content"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        # 输入解析
        if isinstance(data, Path):
            # pathlib.Path → 读取文件内容
            if not data.exists():
                raise FileNotFoundError(f"HTML file not found: {data}")
            self.content = data.read_text(encoding="utf-8")
        elif isinstance(data, str):
            # str → 先检测是否为文件路径，否则作为原始 HTML
            path = Path(data)
            if path.exists() and path.suffix.lower() == ".html":
                self.content = path.read_text(encoding="utf-8")
            else:
                self.content = data
        elif isinstance(data, TransformMedia):
            # Html instances are handled by _unwrap above; this is a type-safety guard
            raise TypeError(
                f"Unsupported HTML data type: {type(data).__name__}. "
                "Please provide a file path (str/Path), raw HTML string, or a file object."
            )
        elif hasattr(data, "read"):
            # TextIO → 读取内容
            data.seek(0)  # type: ignore[union-attr]
            self.content = data.read()  # type: ignore[union-attr]
        else:
            raise TypeError(
                f"Unsupported HTML data type: {type(data).__name__}. "
                "Please provide a file path (str/Path), raw HTML string, or a file object."
            )

        self.caption = caption

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_HTML

    def transform(self, *, step: int, path: Path) -> MediaItem:
        content_encode = self.content.encode("utf-8")
        sha256 = hashlib.sha256(content_encode).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.html"
        fs.safe_write(path / filename, self.content)
        return MediaItem(filename=filename, sha256=sha256, size=len(content_encode), caption=self.caption or "")
