#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 15:54:55
@File: test/test_database_create.py
@IDE: vscode
@Description:
    开发和测试本地数据库的读写能力，并且建立表单
"""
import random
import time
import json

# TODO 在此处导入swanlab中的数据库模块
from swanlab.database.project import SwanProject


# 连接此项目的数据库
sw = SwanProject()

# 迭代次数
epochs = 500
# 学习率
lr = 0.01
# 随机偏移量
offset = random.random() / 5

# FIXME 创建一个实验
sw.init(
    description="this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)

# 模拟训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    # FIXME 在此处将数据写入数据库
    sw.add(tag="loss", data=loss)
    # FIXME 在此处将数据写入数据库
    sw.add(tag="accuracy", data=acc, namespace="train")
    time.sleep(0.1)
