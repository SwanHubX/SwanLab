"""
@author: cunyue
@file: allow.py
@time: 2025/7/3 14:11
@description: 测试 allow 模式的 resume 功能
"""

import random

import swanlab

# 1. 随机生成 run_id
run_id = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=21))


# 2. 初始化第一个实验
# 开启实验并记录
run = swanlab.init(id=run_id, resume='allow')
assert run.id == run_id, "Run ID does not match the expected value"
swanlab.log({"loss": 0.1, "accuracy": 0.9}, step=1)
import time

time.sleep(5)
# 3. 继续实验
run = swanlab.init(id=run_id, resume='allow', reinit=True)
# 上传相同 step 的指标，此时报错
ll = run.log({"loss": 0.1, "accuracy": 0.9}, step=1)
assert ll["loss"].is_error, "Expected loss metric to be in error state due to duplicate step"
assert ll['accuracy'].is_error is not None, "Expected column error to be present for duplicate step"
ll = run.log({"loss": 0.2, "accuracy": 0.8})
assert not ll["loss"].is_error, "Expected loss metric to be logged successfully after reinit"
assert not ll['accuracy'].is_error, "Expected accuracy metric to be logged successfully after reinit"
assert ll['loss'].data == 0.2
assert ll['accuracy'].data == 0.8
