"""
@author: cunyue
@file: test_text.py
@time: 2026/3/11 20:11
@description: 文本处理模块单元测试
"""

import hashlib
from pathlib import Path

from swanlab.proto.swanlab.data.v1.text_pb2 import TextItem, TextValue

# 请根据你的实际路径调整导入
from swanlab.sdk.internal.run.data.transforms.text import Text


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
        """测试 transform 方法的文件名生成与端到端写入逻辑"""
        key = "train/loss_text"
        step = 42
        content = "Sample text for testing"
        caption = "Test caption"

        # 1. 预期计算结果
        expected_sha256 = hashlib.sha256(content.encode()).hexdigest()[:8]
        expected_filename = f"{key}-{step:03d}-{expected_sha256}.__swanlab__.txt"

        # 2. 执行 transform
        result = Text.transform(key=key, step=step, path=tmp_path, content=content, caption=caption)

        # 3. 校验返回的 Protobuf 结构
        assert isinstance(result, TextValue)
        assert len(result.items) == 1

        pb_item: TextItem = result.items[0]
        assert pb_item.filename == expected_filename
        assert pb_item.caption == caption

        # 4. 校验真实文件落盘 (端到端断言)
        target_file = tmp_path / expected_filename
        assert target_file.exists(), f"预期的文件 {target_file} 未在磁盘上生成！"
        assert target_file.read_text(encoding="utf-8") == content

    def test_text_transform_with_nested_input(self, tmp_path: Path):
        """测试 transform 直接接收嵌套的 Text 对象"""
        key = "eval/prompt"
        step = 1
        inner_content = "AI generated text"

        # 构造嵌套输入
        inner_text = Text(content=inner_content, caption="old caption")

        # 传入 transform 时，使用新的 caption 覆盖
        result = Text.transform(key=key, step=step, path=tmp_path, content=inner_text, caption="new overriding caption")

        pb_item: TextItem = result.items[0]
        # 验证套娃解包在 transform 中也生效了
        assert pb_item.caption == "new overriding caption"

        # 验证哈希和文件名是基于解包后的 content 生成的
        expected_sha256 = hashlib.sha256(inner_content.encode()).hexdigest()[:8]
        expected_filename = f"{key}-001-{expected_sha256}.__swanlab__.txt"
        assert pb_item.filename == expected_filename

        # 验证真实文件落盘
        target_file = tmp_path / expected_filename
        assert target_file.exists(), "嵌套提取后的内容未能正确落盘！"
        assert target_file.read_text(encoding="utf-8") == inner_content
