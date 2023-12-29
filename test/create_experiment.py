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


class Enlarge1000Type(sw.BaseType):
    def get_data(self):
        return int(self.value * 1000)

    def get_chart_type(self) -> str:
        return self.chart.line

    def get_namespace(self) -> str:
        return "custom"

    def get_config(self) -> dict:
        return {"color": "#000000"}


# 迭代次数
epochs = 5
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
    sw.log({"test/loss": loss, "accuracy": acc})
    sw.log({"loss2": loss, "accuracy2": acc}, step=epochs - epoch)
    sw.log({"train/loss3": loss, "accuracy3": acc}, step=epoch * 2)
    sw.log({"int": Enlarge1000Type(loss)}, step=epoch * 3)
    time.sleep(0.5)
    # if epoch % 40 == 0:
    #     epoch / 0

print("")
print("")
print("finish training")
