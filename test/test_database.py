#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 15:54:55
@File: test/test_database.py
@IDE: vscode
@Description:
    开发和测试本地数据库的读写能力
"""
import random

# TODO 在此处导入swanlab中的数据库模块


# 一百万次迭代
epochs = 100000
# 学习率
lr = 0.01
# 随机偏移量
offset = random.random() / 5
print(f"lr: {lr}")

# 模拟训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    data = {"accuracy": acc, "loss": loss}
    # TODO 在此处将数据写入数据库
