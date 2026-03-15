"""
@author: cunyue
@file: test_video.py
@time: 2026/3/15
@description: 视频处理模块单元测试
"""

import hashlib
from io import BytesIO
from pathlib import Path

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.media.video_pb2 import VideoItem
from swanlab.sdk.internal.run.transforms.video import Video

# 最小合法 GIF89a（1×1 像素）
GIF_BYTES = (
    b"GIF89a\x01\x00\x01\x00\x80\x00\x00\xff\xff\xff\x00\x00\x00"
    b"!\xf9\x04\x00\x00\x00\x00\x00,\x00\x00\x00\x00\x01\x00\x01\x00\x00\x02\x02D\x01\x00;"
)


@pytest.fixture
def gif_file(tmp_path: Path) -> str:
    """写一个临时 GIF 文件，返回其路径字符串"""
    path = tmp_path / "test.gif"
    path.write_bytes(GIF_BYTES)
    return str(path)


# ---------------------------------- 构造测试 ----------------------------------


class TestVideoInit:
    def test_from_bytes(self):
        v = Video(GIF_BYTES)
        assert v.format == "gif"
        assert v.buffer.getvalue() == GIF_BYTES
        assert v.caption is None

    def test_from_bytesio(self):
        v = Video(BytesIO(GIF_BYTES))
        assert v.format == "gif"
        assert v.buffer.getvalue() == GIF_BYTES

    def test_from_path(self, gif_file):
        v = Video(gif_file)
        assert v.format == "gif"
        assert v.buffer.getvalue() == GIF_BYTES

    def test_caption_stored(self):
        v = Video(GIF_BYTES, caption="my clip")
        assert v.caption == "my clip"

    def test_caption_none_by_default(self):
        v = Video(GIF_BYTES)
        assert v.caption is None


# ---------------------------------- 错误处理测试 ----------------------------------


class TestVideoInitErrors:
    def test_invalid_extension_raises(self, tmp_path):
        p = tmp_path / "clip.mp4"
        p.write_bytes(b"\x00" * 16)
        with pytest.raises(TypeError, match="Unsupported file extension"):
            Video(str(p))

    def test_nonexistent_path_raises(self):
        with pytest.raises(ValueError, match="Failed to open file"):
            Video("/nonexistent/path/clip.gif")

    def test_bad_magic_in_file_raises(self, tmp_path):
        """文件扩展名合法但内容不是 GIF"""
        p = tmp_path / "fake.gif"
        p.write_bytes(b"\xff\xd8\xff\xe0" + b"\x00" * 16)  # JPEG 魔数
        with pytest.raises(TypeError, match="magic number"):
            Video(str(p))

    def test_bad_magic_in_bytes_raises(self):
        with pytest.raises(TypeError, match="Cannot detect video format"):
            Video(b"\x00\x01\x02\x03\x04\x05\x06\x07")

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError, match="Unsupported type"):
            Video(12345)  # type: ignore


# ---------------------------------- 套娃加载测试 ----------------------------------


class TestVideoNesting:
    def test_wrap_copies_buffer_and_format(self):
        inner = Video(GIF_BYTES, caption="inner")
        outer = Video(inner)
        assert outer.buffer is inner.buffer
        assert outer.format == inner.format
        assert outer.caption == "inner"

    def test_outer_caption_overrides_inner(self):
        inner = Video(GIF_BYTES, caption="inner")
        outer = Video(inner, caption="outer")
        assert outer.caption == "outer"

    def test_inner_caption_used_when_outer_none(self):
        inner = Video(GIF_BYTES, caption="inner")
        outer = Video(inner, caption=None)
        assert outer.caption == "inner"


# ---------------------------------- column_type / build_data_record 测试 ----------------------------------


class TestVideoColumnType:
    def test_column_type(self):
        from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType

        assert Video.column_type() == ColumnType.COLUMN_TYPE_VIDEO


class TestVideoBuildDataRecord:
    def test_build_data_record_structure(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Video.build_data_record(key="rollout", step=1, timestamp=ts, data=[item])

        assert record.key == "rollout"
        assert record.step == 1
        assert len(record.videos.items) == 1
        assert record.videos.items[0].filename == item.filename

    def test_build_data_record_multiple_items(self, tmp_path):
        i1 = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        i2 = Video(BytesIO(GIF_BYTES)).transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Video.build_data_record(key="k", step=1, timestamp=ts, data=[i1, i2])
        assert len(record.videos.items) == 2


# ---------------------------------- transform 特有字段测试 ----------------------------------


class TestVideoTransform:
    def test_transform_returns_video_item(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        assert isinstance(item, VideoItem)

    def test_transform_sha256_correct(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        assert item.sha256 == hashlib.sha256(content).hexdigest()

    def test_transform_size_correct(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        assert item.size == len((tmp_path / item.filename).read_bytes())

    def test_transform_caption_empty_when_none(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=1, path=tmp_path)
        assert item.caption == ""

    def test_transform_caption_preserved(self, tmp_path):
        item = Video(GIF_BYTES, caption="hello").transform(step=1, path=tmp_path)
        assert item.caption == "hello"

    def test_transform_format_in_filename(self, tmp_path):
        item = Video(GIF_BYTES).transform(step=3, path=tmp_path)
        assert item.filename.endswith(".gif")
