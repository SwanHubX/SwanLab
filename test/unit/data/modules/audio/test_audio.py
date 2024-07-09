#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:43
@File: pytest_audio.py
@IDE: pycharm
@Description:
    测试音频模块
"""
import os.path

from swanlab.data.modules import Audio
from tutils import TEMP_PATH
from nanoid import generate
import numpy as np
import soundfile as sf
import pytest


def test_audio_ok():
    # ---------------------------------- numpy输入 ----------------------------------

    np.random.randn(2, 100000)
    mock = np.random.randn(2, 100000)
    audio = Audio(data_or_path=mock, sample_rate=44100)
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".wav")
    assert buffer is not None
    assert audio.get_more() is None
    assert audio.get_config() is None

    # ---------------------------------- 路径输入 ----------------------------------

    mock = np.random.randn(2, 100000)
    path = os.path.join(TEMP_PATH, f"{generate()}.wav")
    # 写入文件
    sample_rate = 44200
    sf.write(path, mock.T, sample_rate)
    audio = Audio(data_or_path=path, sample_rate=sample_rate)
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".wav")
    assert buffer is not None


def test_audio_caption():
    # ---------------------------------- numpy输入 ----------------------------------

    np.random.randn(2, 100000)
    mock = np.random.randn(2, 100000)
    audio = Audio(data_or_path=mock, sample_rate=44100, caption="test")
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".wav")
    assert buffer is not None
    assert audio.get_more()["caption"] == "test"
    assert audio.get_config() is None

    # ---------------------------------- 路径输入 ----------------------------------

    mock = np.random.randn(2, 100000)
    path = os.path.join(TEMP_PATH, f"{generate()}.wav")
    # 写入文件
    sample_rate = 44200
    sf.write(path, mock.T, sample_rate)
    audio = Audio(data_or_path=path, sample_rate=sample_rate, caption="test")
    data, buffer = audio.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".wav")
    assert buffer is not None
    assert audio.get_more()["caption"] == "test"
    assert audio.get_config() is None


@pytest.mark.parametrize("dtype", Audio.SF_SUPPORT_DTYPE)
def test_audio_numpy(dtype):
    """
    测试不同的numpy输入情况
    """
    mock = np.random.randn(2, 100000).astype(dtype)
    Audio(data_or_path=mock, sample_rate=44100)
    # 单通道
    mock = np.random.randn(100000).astype(dtype)
    Audio(data_or_path=mock, sample_rate=44100)


def test_audio_fail():
    # 错误的路径
    with pytest.raises(ValueError):
        Audio(data_or_path="not_exist.wav", sample_rate=44100)
    # 不是音频
    path = os.path.join(TEMP_PATH, f"{generate()}.wav")
    with open(path, "w") as f:
        f.write("hello")
    with pytest.raises(ValueError):
        Audio(data_or_path=path, sample_rate=44100)

    # ---------------------------------- 矩阵错误 ----------------------------------

    # 传入矩阵时没有传入采样率
    mock = np.random.randn(2, 100000)
    with pytest.raises(TypeError):
        Audio(data_or_path=mock, sample_rate=None)  # noqa

    # 错误的矩阵类型
    mock = np.random.randn(2, 100000).astype(np.int8)
    with pytest.raises(TypeError):
        Audio(data_or_path=mock, sample_rate=44100)

    # 错误的矩阵形状
    mock = np.random.randn(3, 100000)
    with pytest.raises(TypeError):
        Audio(data_or_path=mock, sample_rate=44100)
