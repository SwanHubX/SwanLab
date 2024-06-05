#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:43
@File: pytest_image.py
@IDE: pycharm
@Description:
    测试图片模块
"""
import os.path

from swanlab.data.modules import Image
from tutils import TEMP_PATH
from nanoid import generate
import numpy as np
import pytest
import torch
from PIL import Image as PILImage
from matplotlib import pyplot as plt


def test_image_ok():
    # ---------------------------------- numpy输入 ----------------------------------

    np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
    mock = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
    image = Image(data_or_path=mock)
    data, buffer = image.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more() is None
    assert image.get_config() is None

    # ---------------------------------- pil输入 ----------------------------------

    mock = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    image = Image(data_or_path=mock)
    data, buffer = image.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more() is None
    assert image.get_config() is None

    # ---------------------------------- 路径输入 ----------------------------------
    mock = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.png")
    # 写入文件
    mock.save(path)
    audio = Image(data_or_path=path)
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None

    # ---------------------------------- pytorch tensor输入 ----------------------------------
    mock = torch.randn(4, 3, 256, 256)
    image = Image(data_or_path=mock)
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more() is None
    assert image.get_config() is None

    # ---------------------------------- plt输入 ----------------------------------
    x = [1, 2, 3]
    y = [2, 3, 5]
    plt.plot(x, y)
    image = Image(data_or_path=plt)
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more() is None
    assert image.get_config() is None


def test_image_caption():
    # ---------------------------------- numpy输入 ----------------------------------

    np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
    mock = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
    image = Image(data_or_path=mock, caption="test")
    data, buffer = image.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more()["caption"] == "test"
    assert image.get_config() is None

    # ---------------------------------- pil输入 ----------------------------------

    mock = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    image = Image(data_or_path=mock, caption="test")
    data, buffer = image.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more()["caption"] == "test"
    assert image.get_config() is None

    # ---------------------------------- 路径输入 ----------------------------------
    mock = PILImage.fromarray(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8))
    path = os.path.join(TEMP_PATH, f"{generate()}.png")
    # 写入文件
    mock.save(path)
    audio = Image(data_or_path=path, caption="test")
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more()["caption"] == "test"
    assert image.get_config() is None

    # ---------------------------------- pytorch tensor输入 ----------------------------------
    mock = torch.randn(4, 3, 256, 256)
    image = Image(data_or_path=mock, caption="test")
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more()["caption"] == "test"
    assert image.get_config() is None

    # ---------------------------------- plt输入 ----------------------------------
    x = [1, 2, 3]
    y = [2, 3, 5]
    plt.plot(x, y)
    image = Image(data_or_path=plt, caption="test")
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".png")
    assert buffer is not None
    assert image.get_more()["caption"] == "test"
    assert image.get_config() is None


def test_image_fail():
    # 错误的路径
    with pytest.raises(ValueError):
        Image(data_or_path="not_exist.png")
    # 不是图片
    path = os.path.join(TEMP_PATH, "test.png")
    with open(path, "w") as f:
        f.write("hello")
    with pytest.raises(ValueError):
        Image(data_or_path=path)

    # ---------------------------------- 矩阵错误 ----------------------------------

    # 通道数错误
    mock = np.random.randint(low=0, high=256, size=(100, 100, 2), dtype=np.uint8)
    with pytest.raises(TypeError):
        Image(data_or_path=mock)

    mock = np.random.randint(low=0, high=256, size=(100, 100, 5), dtype=np.uint8)
    with pytest.raises(TypeError):
        Image(data_or_path=mock)

    # 错误的矩阵形状
    mock = np.random.randint(low=0, high=256, size=(3, 100, 100), dtype=np.uint8)
    with pytest.raises(TypeError):
        Image(data_or_path=mock)

    mock = np.random.randint(low=0, high=256, size=(1, 100, 100, 3), dtype=np.uint8)
    with pytest.raises(TypeError):
        Image(data_or_path=mock)

    # ---------------------------------- 文件类型错误 ----------------------------------
    mock = np.random.randint(low=0, high=256, size=(100, 100, 2), dtype=np.uint8)
    with pytest.raises(ValueError):
        Image(mock, file_type="svg")

    # ---------------------------------- size错误 -------------------------------------
    mock = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
    with pytest.raises(ValueError):
        Image(mock, size=(100, 100, 100))
    with pytest.raises(ValueError):
        Image(mock, size="hello")  # noqa


@pytest.mark.parametrize("file_type", Image.ACCEPT_FORMAT)
def test_image_file_type(file_type):
    """
    测试不同的file_type输入情况
    """
    mock = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)

    image = Image(mock, file_type=file_type)
    assert image.format == file_type

    data, buffer = image.parse()

    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(f".{file_type}")
    assert buffer is not None
    assert image.get_more() is None
    assert image.get_config() is None


def test_image_size():
    """
    测试不同的size输入情况
    """
    mock = np.random.randint(low=0, high=256, size=(256, 512, 3), dtype=np.uint8)
    # 转为PIL图像后，size为(512, 256)

    image = Image(mock, size=None)
    assert image.image_size == (512, 256)

    image = Image(mock, size=(128, 128))
    assert image.image_size == (128, 128)

    image = Image(mock, size=(128, None))
    assert image.image_size == (128, 64)

    image = Image(mock, size=(None, 128))
    assert image.image_size == (256, 128)

    image = Image(mock, size=128)
    assert image.image_size == (128, 64)
