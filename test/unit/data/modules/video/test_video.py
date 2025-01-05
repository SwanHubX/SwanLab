#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/1/5 15:43
@File: pytest_video.py
@IDE: pycharm
@Description:
    测试视频模块
"""
import os.path
import os
import numpy as np
import pytest
import torch
from nanoid import generate

from swanlab.data.modules import Video
from tutils import TEMP_PATH


def test_video_ok():
    # ---------------------------------- numpy输入 ----------------------------------
    
    mock = np.random.randint(low=0, high=256, size=(10, 3, 100, 100), dtype=np.uint8)
    video = Video(data_or_path=mock)
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".gif")
    assert buffer is not None
    assert video.get_more() is None

    # ---------------------------------- 路径输入(带格式) ----------------------------------
    mock = np.random.randint(low=0, high=256, size=(10, 3, 100, 100), dtype=np.uint8)
    path = os.path.join(TEMP_PATH, f"{generate()}.mp4")
    # 写入文件
    with open(path, "wb") as f:
        f.write(mock)
    video = Video(data_or_path=path, format="mp4")
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".mp4")
    assert buffer is not None
    assert video.get_more() is None
    
    # ---------------------------------- 路径输入(不带格式) ----------------------------------
    mock = np.random.randint(low=0, high=256, size=(10, 3, 100, 100), dtype=np.uint8)
    path = os.path.join(TEMP_PATH, f"{generate()}.mp4")
    # 写入文件
    with open(path, "wb") as f:
        f.write(mock)
    video = Video(data_or_path=path)
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".gif")
    assert buffer is not None
    assert video.get_more() is None

    # ---------------------------------- pytorch tensor输入 ----------------------------------
    mock = torch.randint(
        low=0, high=256, size=(10, 3, 100, 100), dtype=torch.uint8
    )
    video = Video(data_or_path=mock, format="mp4")
    data, buffer = video.parse()
    # 返回文件名称
    assert isinstance(data, str)
    assert data.endswith(".mp4")
    assert buffer is not None
    assert video.get_more() is None