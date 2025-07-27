"""
@author: cunyue
@file: must.py
@time: 2025/7/3 14:11
@description: 测试 must 模式的 resume 功能
"""

import random

import swanlab

# 1. 初始化第一个实验
# 开启实验并记录
run = swanlab.init()
swanlab.log({"loss": 0.1, "accuracy": 0.9}, step=1)

import time

time.sleep(5)
# 3. 继续第一个实验
run = swanlab.init(id=run.id, resume='must', reinit=True)
# 上传相同 step 的指标，此时报错
ll = run.log({"loss": 0.3, "accuracy": 0.7}, step=1)
assert ll["loss"].is_error, "Expected loss metric to be in error state due to duplicate step"
assert ll['accuracy'].is_error is not None, "Expected column error to be present for duplicate step"
ll = run.log({"loss": 0.5, "accuracy": 0.8})
assert not ll["loss"].is_error, "Expected loss metric to be logged successfully after reinit"
assert not ll['accuracy'].is_error, "Expected accuracy metric to be logged successfully after reinit"
assert ll['loss'].data == 0.5
assert ll['accuracy'].data == 0.8

# 2. 继续一个不存在的实验
try:
    run = swanlab.init(
        id="".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=21)),
        resume='must',
        reinit=True,
    )
    raise Exception("Should have raised exception")
except RuntimeError:
    pass
