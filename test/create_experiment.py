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


class Enlarge100000Billion(sw.data.BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        return self.value * 100000000000000000

    def get_config(self, *args, **kwargs) -> dict:
        return {"color": "red"}

    def get_namespace(self, *args, **kwargs) -> str:
        return "custom"

    def get_chart_type(self) -> str:
        return self.chart.line


class Shrink100Billion(sw.data.BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        return self.value / 100000000000000

    def get_config(self, *args, **kwargs) -> dict:
        return {"color": "red"}

    def get_namespace(self, *args, **kwargs) -> str:
        return "custom"

    def get_chart_type(self) -> str:
        return self.chart.line


# 迭代次数
epochs = 5000
# 学习率
lr = 0.01
# 随机偏移量
offset = random.random() / 5
# 创建一个实验
run = sw.init(
    description=" this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
    logdir="swanlog",
    log_level="debug",
)


print(run.config.learning_rate)

print("start training")

print("")

# 模拟训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    sw.log({"loss": loss, "accuracy": acc}, step=epoch)
    if epoch < 10:
        sw.log(
            {
                "loss_enlarge100000Billion": Enlarge100000Billion(loss),
                "accuracy_shrink100Billion": Shrink100Billion(acc),
            },
            step=1,
        )
    else:
        sw.log(
            {
                "loss_enlarge100000Billion": Enlarge100000Billion(loss),
                "accuracy_shrink100Billion": Shrink100Billion(acc),
            },
            step=epoch,
        )
    # sw.log({"accuracy2": f"{acc}", "test/loss2": f"is {loss}"}, step=epochs - epoch)
    # sw.log({"loss3": loss, "accuracy3": acc}, step=1)
    # sw.log({"loss4": loss, "accuracy4": acc}, step=epoch * 2)
    time.sleep(0.05)

print("")
print("")
print("finish training")
