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
import numpy as np
from PIL import Image as PILImage
import os


class Enlarge1000(sw.data.BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        return int(self.value * 1000)

    def get_config(self, *args, **kwargs) -> dict:
        return {"color": "red"}

    def get_namespace(self, *args, **kwargs) -> str:
        return "custom"

    def get_chart_type(self) -> str:
        return self.chart.image


# 迭代次数
epochs = 3
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

    current_file_path = os.path.dirname(os.path.abspath(__file__))
    print("current_file_path:", current_file_path)
    test_path = os.path.join(current_file_path, "assets/test_image.png")

    test_pil = PILImage.new("RGB", (100, 100), color="red")
    test_nparray = np.random.randint(255, size=(300, 300, 3), dtype=np.uint8)

    # if epoch < 10:
    #     sw.log({"loss": Enlarge1000(loss), "accuracy": acc}, step=1)
    # else:
    # sw.log({"loss": sw.data.Image("./test_image.jpg"), "accuracy": acc})

    sw.log(
        {
            "String": sw.Image(test_path),
            "PILImage": sw.Image(test_pil),
            "nparray": sw.Image(test_nparray),
            "accuracy": acc,
        }
    )
    # sw.log({"loss": Enlarge1000(loss), "accuracy": acc})

    # sw.log({"accuracy2": f"{acc}", "test/loss2": f"is {loss}"}, step=epochs - epoch)
    # sw.log({"loss3": loss, "accuracy3": acc}, step=1)
    # sw.log({"loss4": loss, "accuracy4": acc}, step=epoch * 2)
    time.sleep(0.5)

print("")
print("")
print("finish training")
