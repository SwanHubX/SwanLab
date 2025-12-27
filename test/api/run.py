"""
@author: Zhou QiYang
@file: run.py
@time: 2025/12/27 15:44
@description: 获取单个实验功能的单元测试
"""

import wandb

import swanlab


def test_example_code():
    api = swanlab.OpenApi()
    exp = api.run(path="bainiantest/SwanLab122/o6b0auazzj9vxxyu6rjd1")  # 可通过api.runs()获取expid
    print(exp.history(keys=['loss'], sample=20, x_axis='t/accuracy', pandas=True))

    wandb.Api
