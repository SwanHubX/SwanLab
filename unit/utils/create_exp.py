#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-01 21:22:19
@File: test/unit/utils/create_exp.py
@IDE: vscode
@Description:
    创建实验，会自动创建unit/temp目录，存放实验数据
"""
import os

# 当前文件的绝对路径
cur_path = os.path.abspath(__file__)
# swanlog存储位置
save_dir = os.path.join(os.path.dirname(os.path.dirname(cur_path)), "temp")
# swanlab库存储位置
swanlab_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(cur_path))))
import sys

sys.path.append(swanlab_dir)
import numpy as np
from swanlab.env import ROOT
import swanlab

os.environ[ROOT] = save_dir
print("set swanlog path: ", os.environ[ROOT])


if not os.path.exists(save_dir):
    os.mkdir(save_dir)
run = swanlab.init(log_level="debug")
# audio
sample_rate = 44100
test_audio_arr = np.random.randn(2, 100000)
run.log({"test-1": 1, "test-2": 1, "no": swanlab.Audio(test_audio_arr, sample_rate=sample_rate)})
swanlab.finish()


if __name__ == "__main__":
    # 注册的实验名称
    exp_name = run.settings.exp_name
    # 检查写入的字段是否正确,不需要再次连接数据库
    from swanlab.db import *

    # 数据库中应该只有一个名为exp_name的实验
    assert Experiment.select().where(Experiment.name == exp_name).count() == 1
    # 获取这个实验
    exp = Experiment.get(name=exp_name)
    assert exp is not None
    assert exp.run_id == run.settings.run_id
    assert exp.name == exp_name
    assert exp.status == 1
    # 新建实验的sort应该等于实验的id（目前是这样的
    assert exp.sort == exp.id
    """
    检查当前实验的tag字段数量
    """
    print(exp.tags)

    print(exp)
