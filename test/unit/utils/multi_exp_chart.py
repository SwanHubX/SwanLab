#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-01 17:35:21
@File: test/unit/utils/multi_exp_chart.py
@IDE: vscode
@Description:
    测试多实验图表后端数据生成
"""
from swanlab.env import ROOT
import swanlab
import os
import numpy as np
import shutil

# 当前文件的绝对路径
cur_path = os.path.abspath(__file__)
# swanlog存储位置
save_dir = os.path.join(os.path.dirname(os.path.dirname(cur_path)), "temp", "mutli_exp_chart")


def run():
    """运行试验"""
    # 设置swanlog存储位置
    os.environ[ROOT] = save_dir
    if os.path.exists(os.environ[ROOT]):
        shutil.rmtree(os.environ[ROOT])
    os.mkdir(os.environ[ROOT])
    print("set swanlog path: ", os.environ[ROOT])
    run = swanlab.init(log_level="debug")
    # audio
    sample_rate = 44100
    test_audio_arr = np.random.randn(2, 100000)
    run.log({"test-1": 1, "test-2": 1, "no": swanlab.Audio(test_audio_arr, sample_rate=sample_rate)})


if __name__ == "__main__":
    # 开启两个进程，写入两个实验
    import multiprocessing

    process1 = multiprocessing.Process(target=run)
    process2 = multiprocessing.Process(target=run)

    process1.start()
    process2.start()

    process1.join()
    process2.join()

    # 开始校验数据库信息
    from swanlab.db import *

    # 连接数据库
    connect()
    # ...
