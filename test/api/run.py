"""
@author: Zhou QiYang
@file: run.py
@time: 2025/12/27 15:44
@description: 获取单个实验功能的单元测试
"""

import time

import swanlab


def test_history_example_code():
    api = swanlab.Api()
    exp = api.run(path="username/project/expid")  # 可通过api.runs()获取expid
    start = time.perf_counter()

    print(exp.history())
    end = time.perf_counter()
    print(f"total time: {end - start} seconds")

    print(exp.history(keys=exp.metric_keys))
    end = time.perf_counter()
    print(f"total time: {end - start} seconds")
