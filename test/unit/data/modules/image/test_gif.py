"""
@author: cunyue
@file: test_gif.py
@time: 2025/7/17 11:44
@description: 测试 Image 模块的GIF功能
"""

import pytest

from swanlab.data.modules.image import Image
from swanlab.toolkit import MediaBuffer


class TestImageGif:
    @staticmethod
    def test_create_image_from_gif_path(tmp_path):
        gif_path = tmp_path / "sample.gif"
        gif_path.write_bytes(b"GIF89a")  # 简单的GIF头部
        img = Image(str(gif_path))
        assert img.format == "gif"
        assert isinstance(img.buffer, MediaBuffer)
        assert img.caption is None or isinstance(img.caption, str)

    @staticmethod
    def test_create_image_from_gif_path_with_caption(tmp_path):
        gif_path = tmp_path / "sample.gif"
        gif_path.write_bytes(b"GIF89a")
        img = Image(str(gif_path), caption="hello world")
        assert img.format == "gif"
        assert img.get_more() == {"caption": "hello world"}

    @staticmethod
    def test_invalid_gif_path_raises():
        with pytest.raises(FileNotFoundError):
            Image("not_exist.gif")
