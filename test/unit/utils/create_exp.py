#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-01 21:22:19
@File: test/unit/utils/create_exp.py
@IDE: vscode
@Description:
    创建实验
"""
import os

# 当前文件的绝对路径
cur_path = os.path.abspath(__file__)
# swanlog存储位置
save_dir = os.path.join(os.path.dirname(os.path.dirname(cur_path)), "temp", "test")
# swanlab库存储位置
swanlab_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(cur_path))))
import sys

sys.path.append(swanlab_dir)
import numpy as np
from swanlab.env import ROOT
import swanlab

os.environ[ROOT] = save_dir
print("set swanlog path: ", os.environ[ROOT])


def run():
    """运行实验"""
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    run = swanlab.init(log_level="debug")
    # audio
    sample_rate = 44100
    test_audio_arr = np.random.randn(2, 100000)
    run.log({"test-1": 1, "test-2": 1, "no": swanlab.Audio(test_audio_arr, sample_rate=sample_rate)})
    try:
        swanlab.finish(another_run=True)
    except RuntimeError:
        pass


if __name__ == "__main__":
    run()
