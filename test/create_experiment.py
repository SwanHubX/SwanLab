#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 15:54:55
@File: test/create_experiment.py
@IDE: vscode
@Description:
    开启一个实验
"""
import random
import swanlab as sw
import time

# 迭代次数
epochs = 100
# 学习率
lr = 0.01
# 随机偏移量
offset = random.random() / 5
# 创建一个实验
sw.init(
    description="this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
    log_level="debug",
)

print("start training")

print("")

# 模拟训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    sw.log({"loss": loss, "accuracy": acc})
    sw.log({"loss2": loss, "accuracy2": acc}, step=epochs - epoch)
    sw.log({"loss3": loss, "accuracy3": acc}, step=epoch * 2)

    time.sleep(0.5)
    # if epoch % 40 == 0:
    #     epoch / 0

print("")
print("")
print("finish training")
