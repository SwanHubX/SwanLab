"""
@author: caddiesnew
@file: test_html.py
@time: 2026/6/8
@description: HTML TransformMedia 单元测试
"""

import hashlib
from io import BytesIO, StringIO
from pathlib import Path

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.run.transforms.html import Html

# ---------------------------------- 构造测试 ----------------------------------


class TestHtmlInit:
    def test_from_raw_string(self):
        """从原始 HTML 字符串构造"""
        h = Html("<h1>Hello</h1><p>World</p>", caption="greeting")
        assert h.content == "<h1>Hello</h1><p>World</p>"
        assert h.caption == "greeting"

    def test_from_raw_string_plain_text(self):
        """非 .html 后缀的 str 即使文件不存在也作为原始字符串"""
        h = Html("just some text, not a file")
        assert h.content == "just some text, not a file"

    def test_from_html_file_string(self, tmp_path: Path):
        """str 传入存在的 .html 文件路径 → 读取文件内容"""
        html_file = tmp_path / "report.html"
        html_file.write_text("<html><body>file content</body></html>", encoding="utf-8")
        h = Html(str(html_file))
        assert h.content == "<html><body>file content</body></html>"

    def test_from_nonexistent_html_path_string(self):
        """str 传入不存在的 .html 路径 → 视为原始 HTML 字符串"""
        h = Html("./not_exist.html")
        assert h.content == "./not_exist.html"

    def test_from_path_object(self, tmp_path: Path):
        """pathlib.Path 输入 → 读取文件内容"""
        html_file = tmp_path / "page.html"
        html_file.write_text("<div>path object</div>", encoding="utf-8")
        h = Html(html_file)
        assert h.content == "<div>path object</div>"

    def test_from_path_object_not_found(self, tmp_path: Path):
        """pathlib.Path 指向不存在的文件 → FileNotFoundError"""
        missing = tmp_path / "missing.html"
        with pytest.raises(FileNotFoundError, match="HTML file not found"):
            Html(missing)

    def test_from_file_object(self):
        """文件类对象 (TextIO) 输入"""
        f = StringIO("<p>from file object</p>")
        h = Html(f)
        assert h.content == "<p>from file object</p>"

    def test_from_file_object_with_seek(self):
        """文件对象 seek(0) 后读取完整内容"""
        f = StringIO("<span>seek test</span>")
        f.seek(5)  # 移动指针
        h = Html(f)
        assert h.content == "<span>seek test</span>"

    def test_caption_none_by_default(self):
        """caption 默认为 None"""
        h = Html("<p>no caption</p>")
        assert h.caption is None

    def test_empty_html_string(self):
        """空 HTML 字符串"""
        h = Html("")
        assert h.content == ""

    def test_unicode_html(self):
        """含 Unicode 的 HTML"""
        content = "<h1>你好世界 🌍</h1>"
        h = Html(content)
        assert h.content == content

    def test_long_raw_html_string_not_treated_as_path(self):
        """超长原始 HTML 字符串 (超过路径上限) 即使以 .html 结尾也不触发文件系统调用

        修复前: Path(超长字符串).is_file() 会抛 OSError (File name too long) 崩掉用户 run
        """
        # 长度超过 Linux PATH_MAX (4096), 且刻意以 .html 结尾以触发「路径 vs 内容」判定
        content = "<body>" + "x" * 5000 + "</body>.html"
        h = Html(content)
        assert h.content == content

    def test_html_path_string_with_null_byte(self):
        """以 .html 结尾但含 null byte 的字符串走 try-except 兜底, 作为原始内容

        修复前: Path(含 null).is_file() 抛 ValueError (embedded null byte)
        """
        content = "report\x00.html"
        h = Html(content)
        assert h.content == content

    def test_path_pointing_to_directory_raises_error(self, tmp_path: Path):
        """Path 指向目录 → FileNotFoundError (而非 IsADirectoryError)"""
        d = tmp_path / "subdir"
        d.mkdir()
        with pytest.raises(FileNotFoundError, match="is not a file"):
            Html(d)

    def test_binary_file_object_decoded(self):
        """二进制 file-like (BytesIO) 内容被解码为 str, 避免 transform() 中 encode() 失败"""
        f = BytesIO("<p>binary mode</p>".encode("utf-8"))
        h = Html(f)
        assert h.content == "<p>binary mode</p>"

    def test_non_seekable_file_object(self):
        """不可 seek 的 file-like 对象不崩, 仍能读取内容"""

        class NonSeekableStream:
            def read(self) -> str:
                return "<p>non-seekable</p>"

            def seek(self, *_) -> int:
                raise OSError("stream is not seekable")

        h = Html(NonSeekableStream())  # type: ignore[arg-type]
        assert h.content == "<p>non-seekable</p>"


# ---------------------------------- 套娃加载测试 ----------------------------------


class TestHtmlNesting:
    def test_wrap_copies_content(self):
        """套娃加载复用内层 content"""
        inner = Html("<p>inner</p>", caption="inner caption")
        outer = Html(inner)
        assert outer.content == "<p>inner</p>"
        assert outer.caption == "inner caption"

    def test_outer_caption_overrides_inner(self):
        """外层 caption 优先级高于内层"""
        inner = Html("<p>inner</p>", caption="inner")
        outer = Html(inner, caption="outer")
        assert outer.content == "<p>inner</p>"
        assert outer.caption == "outer"

    def test_inner_caption_used_when_outer_none(self):
        """外层 caption 为 None 时使用内层 caption"""
        inner = Html("<p>inner</p>", caption="inner")
        outer = Html(inner, caption=None)
        assert outer.caption == "inner"


# ---------------------------------- column_type 测试 ----------------------------------


class TestHtmlColumnType:
    def test_column_type(self):
        assert Html.column_type() == ColumnType.COLUMN_TYPE_HTML


# ---------------------------------- transform 测试 ----------------------------------


class TestHtmlTransform:
    def test_transform_returns_media_item(self, tmp_path: Path):
        """transform 返回 MediaItem"""
        item = Html("<h1>Test</h1>").transform(step=1, path=tmp_path)
        assert isinstance(item, MediaItem)

    def test_transform_html_filename(self, tmp_path: Path):
        """输出文件名格式为 {step:03d}-{sha256[:8]}.html"""
        content = "<h1>Hello</h1>"
        item = Html(content).transform(step=5, path=tmp_path)
        expected_sha256 = hashlib.sha256(content.encode("utf-8")).hexdigest()[:8]
        assert item.filename == f"005-{expected_sha256}.html"

    def test_transform_sha256_correct(self, tmp_path: Path):
        """MediaItem.sha256 与内容一致"""
        content = "hello world"
        item = Html(content).transform(step=1, path=tmp_path)
        assert item.sha256 == hashlib.sha256(content.encode("utf-8")).hexdigest()

    def test_transform_size_correct(self, tmp_path: Path):
        """MediaItem.size 与内容字节数一致"""
        content = "<p>size test</p>"
        item = Html(content).transform(step=1, path=tmp_path)
        assert item.size == len(content.encode("utf-8"))

    def test_transform_file_written(self, tmp_path: Path):
        """文件内容正确写入"""
        content = "<html><body>written content</body></html>"
        item = Html(content).transform(step=0, path=tmp_path)
        assert (tmp_path / item.filename).exists()
        assert (tmp_path / item.filename).read_text(encoding="utf-8") == content

    def test_transform_caption_empty_when_none(self, tmp_path: Path):
        """caption 为 None 时返回空字符串"""
        item = Html("<p>no caption</p>").transform(step=1, path=tmp_path)
        assert item.caption == ""

    def test_transform_caption_preserved(self, tmp_path: Path):
        """caption 有值时正确传递"""
        item = Html("<p>with caption</p>", caption="my caption").transform(step=1, path=tmp_path)
        assert item.caption == "my caption"

    def test_transform_unicode_content(self, tmp_path: Path):
        """Unicode 内容正确写入"""
        content = "<h1>你好世界 🌍</h1>"
        item = Html(content).transform(step=0, path=tmp_path)
        assert (tmp_path / item.filename).read_text(encoding="utf-8") == content

    def test_transform_nested_content_written(self, tmp_path: Path):
        """套娃时，落盘内容为内层解包后的 content"""
        inner = Html("<p>nested content</p>", caption="old")
        outer = Html(inner, caption="new")
        result = outer.transform(step=1, path=tmp_path)
        assert result.caption == "new"
        assert (tmp_path / result.filename).read_text(encoding="utf-8") == "<p>nested content</p>"

    def test_transform_file_exists_on_disk(self, tmp_path: Path):
        """生成的 .html 文件确实存在于磁盘"""
        item = Html("<div>exists test</div>").transform(step=10, path=tmp_path)
        assert (tmp_path / item.filename).is_file()


# ---------------------------------- build_data_record 测试 ----------------------------------


class TestHtmlBuildDataRecord:
    def test_build_data_record_structure(self, tmp_path: Path):
        """build_data_record 返回正确结构的 MediaRecord"""
        html = Html("<h1>record test</h1>")
        item = html.transform(step=1, path=tmp_path)
        ts = Timestamp(seconds=1234567890)

        record = Html.build_data_record(key="page", step=1, timestamp=ts, data=[item])

        assert record.key == "page"
        assert record.step == 1
        assert record.timestamp == ts
        assert record.type == ColumnType.COLUMN_TYPE_HTML
        assert len(record.value.items) == 1
        assert record.value.items[0].filename == item.filename
        assert record.value.items[0].sha256 == item.sha256
        assert record.value.items[0].size == item.size

    def test_build_data_record_multiple_items(self):
        """build_data_record 支持多个 MediaItem"""
        items = [
            MediaItem(filename="000-abc12345.html", sha256="hash1", size=100, caption="cap1"),
            MediaItem(filename="001-def67890.html", sha256="hash2", size=200, caption="cap2"),
        ]
        ts = Timestamp()
        record = Html.build_data_record(key="pages", step=0, timestamp=ts, data=items)

        assert record.key == "pages"
        assert len(record.value.items) == 2
        assert record.value.items[0].filename == "000-abc12345.html"
        assert record.value.items[1].caption == "cap2"


# ---------------------------------- 构造错误测试 ----------------------------------


class TestHtmlInitErrors:
    def test_unsupported_type_int(self):
        """传入 int 应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported HTML data type"):
            Html(12345)  # type: ignore

    def test_unsupported_type_list(self):
        """传入 list 应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported HTML data type"):
            Html(["<p>list</p>"])  # type: ignore

    def test_unsupported_type_none(self):
        """传入 None 应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported HTML data type"):
            Html(None)  # type: ignore


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
