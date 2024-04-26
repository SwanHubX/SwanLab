#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-04 13:26:59
@File: test/create_experiment.py
@IDE: vscode
@Description:
    创建一个文件，作为测试用例
    WARNING 请勿随意修改此文件，以免影响测试效果
"""
import swanlab
import time
import random
import numpy as np

epochs = 50
lr = 0.01
offset = random.random() / 5
# 初始化
swanlab.init(
    description="这是一个测试实验",
    log_level="debug",
    config="test/config/config.json",
    load="test/config/load.yaml",
    # cloud=True,
)
swanlab.config.epoches = epochs
swanlab.config.learning_rate = lr
swanlab.config.debug = "这是一串" + "很长" * 100 + "的字符串"
# 模拟训练
for epoch in range(2, swanlab.config.epoches):
    acc = 1 - 2 ** -epoch - random.random() / epoch - offset
    loss = 2 ** -epoch + random.random() / epoch + offset
    loss2 = 3 ** -epoch + random.random() / epoch + offset * 3
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    if epoch % 10 == 0:
        # 测试audio
        sample_rate = 44100
        test_audio_arr = np.random.randn(2, 100000)
        swanlab.log(
            {
                "test/audio": [swanlab.Audio(test_audio_arr, sample_rate, caption="test")] * (epoch // 10),
            },
            step=epoch,
        )
        # 测试image
        test_image = np.random.randint(0, 255, (100, 100, 3))
        swanlab.log(
            {
                "test/image": [swanlab.Image(test_image, caption="test")] * (epoch // 10),
            },
            step=epoch,
        )
        # 测试text
        swanlab.log(
            {
                "text": swanlab.Text("这是一段测试文本", caption="test"),
            },
            step=epoch,
        )
        # 测试折线图
        swanlab.log({"t/accuracy": acc, "loss": loss, "loss2": loss2})
    else:
        # 测试折线图
        swanlab.log({"t/accuracy": acc, "loss": loss, "loss2": loss2})
    time.sleep(0.5)
