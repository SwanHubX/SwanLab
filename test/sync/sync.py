"""
@author: cunyue
@file: sync.py
@time: 2025/7/20 17:23
@description: 测试同步功能，与同目录下的sync.sh文件配合使用，使用方式为：
1. python test/sync/sync.py: 先运行 sync.sh，此时会要求输入 run_dir
2. bash ./test/sync/sync.s: 然后在项目根目录下运行此脚本，他会打印出当前实验的路径，复制此路径到步骤一，按回车
这样将实现每五秒钟同步一次数据到云端
"""

import os
import random
import time

import numpy as np

import swanlab
from tutils import TEMP_PATH

logdir = os.path.join(TEMP_PATH, "sync")

run = swanlab.init(logdir=logdir, mode="offline", project="test-sync", description="test-sync", experiment_name="sync")

input("日志文件夹路径为: " + run.public.run_dir + " , 按任意键进入下一步")

epochs = 500
lr = 0.01
offset = random.random() / 5

swanlab.config.epochs = epochs
swanlab.config.learning_rate = lr
swanlab.config.offset = offset

# 模拟训练
for epoch in range(2, swanlab.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    if epoch % 10 == 0:
        # 测试image
        images = [swanlab.Image(np.random.randint(0, 255, (100, 100, 3)), caption="test") for _ in range(epoch // 10)]
        swanlab.log({"test/image": images}, step=epoch)
        # 测试折线图
        swanlab.log({"t/accuracy": acc, "loss": loss}, step=epoch)
    else:
        # 测试折线图
        swanlab.log({"t/accuracy": acc, "loss": loss}, step=epoch)
    time.sleep(1)
