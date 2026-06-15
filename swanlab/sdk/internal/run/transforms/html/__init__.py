"""
@author: caddiesnew
@file: __init__.py
@time: 2026/6/8
@description: HTML transform module — converts HTML strings or files into protobuf MediaItem records
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

# 视为「可能是文件路径」的字符串长度上限。
# 取 Linux PATH_MAX (4096) —— 主流平台最宽松的合法路径上限; macOS (1024) / Windows (260/32767) 更短的限制
# 由下方 try/except 兜底, 超长路径不会崩溃; 此上限仅用于避免对明显是 HTML 内容的长串做无意义的 stat 调用。
_PATH_LEN_LIMIT = 4096


class Html(TransformMedia):
    def __init__(self, data: HtmlDataType, caption: CaptionType = None):
        """Html class constructor

        Parameters
        ----------
        data: str, pathlib.Path, TextIO, or Html
            - str: an HTML file path (must end with .html and exist on disk) or a raw HTML string.
            - pathlib.Path: path to an HTML file.
            - TextIO: a file-like object, e.g. the return value of open("file.html").
            - Html: nesting — copies attributes from another Html instance.
        caption: str, optional
            Optional caption for the HTML content.
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
            if not data.is_file():
                raise FileNotFoundError(f"HTML file not found or is not a file: {data}")
            self.content = data.read_text(encoding="utf-8")
        elif isinstance(data, str):
            # str → 先检测是否为文件路径，否则作为原始 HTML
            # 注意: 原始 HTML 字符串可能很长或含特殊字符 (如 null byte),
            # 直接 Path(data).exists() 会抛 OSError (File name too long) 或 ValueError, 故先做廉价过滤再兜底
            try:
                if len(data) < _PATH_LEN_LIMIT and data.lower().endswith(".html"):
                    path = Path(data)
                    if path.is_file():
                        self.content = path.read_text(encoding="utf-8")
                    else:
                        self.content = data
                else:
                    self.content = data
            except (OSError, ValueError):
                self.content = data
        elif isinstance(data, TransformMedia):
            # Html instances are handled by _unwrap above; this is a type-safety guard
            raise TypeError(
                f"Unsupported HTML data type: {type(data).__name__}. "
                "Please provide a file path (str/Path), raw HTML string, or a file object."
            )
        elif hasattr(data, "read"):
            # TextIO → 读取内容
            # 不可 seek 的流 (如 socket/管道) 调用 seek(0) 会抛异常, 忽略即可
            try:
                data.seek(0)  # type: ignore[union-attr]
            except (OSError, ValueError):
                pass
            raw = data.read()  # type: ignore[union-attr]
            # 二进制 file-like 的 read() 返回 bytes, 统一解码为 str 避免 transform() 中 encode() 失败
            self.content = raw.decode("utf-8") if isinstance(raw, bytes) else raw
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
