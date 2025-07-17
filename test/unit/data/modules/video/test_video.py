#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:43
@File: test_video.py
@IDE: pycharm
@Description:
    测试视频模块
"""
import os.path

import numpy as np
import pytest
from PIL import Image as PILImage
from nanoid import generate

from swanlab.data.modules import Video
from tutils import TEMP_PATH


def test_video_ok():
    """测试Video模块的基本功能"""
    # ---------------------------------- GIF文件路径输入 ----------------------------------
    # 创建一个测试用的GIF文件
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.gif")
    # 保存为GIF格式
    mock_image.save(path, format="GIF")

    video = Video(data_or_path=path)
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".gif")
    assert buffer is not None
    assert video.get_more() is None


def test_video_caption():
    """测试Video模块的caption功能"""
    # ---------------------------------- GIF文件路径输入带caption ----------------------------------
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.gif")
    mock_image.save(path, format="GIF")

    video = Video(data_or_path=path, caption="test video")
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".gif")
    assert buffer is not None
    assert video.get_more()["caption"] == "test video"


def test_video_fail():
    """测试Video模块的错误处理"""
    # 错误的路径
    with pytest.raises(FileNotFoundError):
        Video(data_or_path="not_exist.gif")

    # 不是GIF格式的文件
    path = os.path.join(TEMP_PATH, "test.txt")
    with open(path, "w") as f:
        f.write("hello")
    with pytest.raises(ValueError):
        Video(data_or_path=path)

    # 不是GIF格式的图片文件
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.png")
    mock_image.save(path)
    with pytest.raises(ValueError):
        Video(data_or_path=path)

    # 其他格式的视频文件（不实际创建文件，只测试文件名验证）
    path = os.path.join(TEMP_PATH, f"{generate()}.mp4")
    with pytest.raises(ValueError):
        Video(data_or_path=path)


def test_video_format_validation():
    """测试Video模块的格式验证"""
    # 测试各种非GIF格式
    non_gif_formats = [
        ".png",
        ".jpg",
        ".jpeg",
        ".bmp",
        ".tiff",
        ".mp4",
        ".avi",
        ".mov",
        ".wmv",
        ".flv",
        ".webm",
        ".mkv",
    ]

    for format_ext in non_gif_formats:
        path = os.path.join(TEMP_PATH, f"{generate()}{format_ext}")
        with pytest.raises(ValueError, match="swanlab.Video only supports gif format file paths"):
            Video(data_or_path=path)


def test_video_parse_consistency():
    """测试Video模块parse方法的一致性"""
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.gif")
    mock_image.save(path, format="GIF")

    video = Video(data_or_path=path, caption="test")

    # 多次调用parse应该返回相同结果
    data1, buffer1 = video.parse()
    data2, buffer2 = video.parse()

    assert data1 == data2
    assert buffer1 == buffer2
    assert data1.endswith(".gif")


def test_video_file_cleanup():
    """测试Video模块的文件清理"""
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.gif")
    mock_image.save(path, format="GIF")

    video = Video(data_or_path=path)
    data, buffer = video.parse()

    # 验证文件存在
    assert os.path.exists(path)
    assert isinstance(data, str)
    assert data.endswith(".gif")
    assert buffer is not None


def test_video_inheritance():
    """测试Video模块是否正确继承了媒体处理的功能"""
    mock_image = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.gif")
    mock_image.save(path, format="GIF")

    video = Video(data_or_path=path, caption="test")
    # 验证Video有Image的所有必要方法
    assert hasattr(video, 'parse')
    assert hasattr(video, 'get_more')
