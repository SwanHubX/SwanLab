"""
@author: cunyue
@file: test_text.py
@time: 2026/3/11 20:11
@description: 文本处理模块单元测试
"""

from pathlib import Path

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.media.text_pb2 import TextItem
from swanlab.sdk.internal.run.transforms.text import Text


class TestTextTransform:
    def test_text_init_basic(self):
        """测试基础的实例化"""
        t = Text(content="hello world", caption="a greeting")
        assert t.content == "hello world"
        assert t.caption == "a greeting"

    def test_text_init_nested_override(self):
        """测试套娃加载：外层参数覆盖内层参数"""
        inner = Text(content="inner content", caption="inner caption")
        # 外层提供了新的 caption，应该覆盖内层的
        outer = Text(content=inner, caption="outer caption")

        assert outer.content == "inner content"
        assert outer.caption == "outer caption"

    def test_text_init_nested_fallback(self):
        """测试套娃加载：外层参数为空时回退使用内层参数"""
        inner = Text(content="inner content", caption="inner caption")
        # 外层未提供 caption (None)，应该保留内层的
        outer = Text(content=inner)

        assert outer.content == "inner content"
        assert outer.caption == "inner caption"

    def test_text_transform(self, tmp_path: Path):
        """transform 将内容正确写入文件，caption 写入 TextItem"""
        text = Text(content="Sample text for testing", caption="Test caption")
        result = text.transform(step=42, path=tmp_path)

        assert isinstance(result, TextItem)
        assert result.caption == "Test caption"
        assert (tmp_path / result.filename).read_text(encoding="utf-8") == "Sample text for testing"

    def test_text_transform_nested_content_written(self, tmp_path: Path):
        """套娃时，落盘内容为内层解包后的 content"""
        inner = Text(content="AI generated text", caption="old caption")
        outer = Text(content=inner, caption="new caption")
        result = outer.transform(step=1, path=tmp_path)

        assert result.caption == "new caption"
        assert (tmp_path / result.filename).read_text(encoding="utf-8") == "AI generated text"

    def test_text_build_data_record(self):
        """测试 build_data_record 方法"""
        key = "test/metric"
        step = 10
        timestamp = Timestamp(seconds=1234567890)

        items = [
            TextItem(filename="file1.txt", caption="caption1"),
            TextItem(filename="file2.txt", caption="caption2"),
        ]

        record = Text.build_data_record(key=key, step=step, timestamp=timestamp, data=items)

        assert record.key == key
        assert record.step == step
        assert record.timestamp == timestamp
        assert len(record.texts.items) == 2
        assert record.texts.items[0].filename == "file1.txt"
        assert record.texts.items[0].caption == "caption1"
        assert record.texts.items[1].filename == "file2.txt"
        assert record.texts.items[1].caption == "caption2"
