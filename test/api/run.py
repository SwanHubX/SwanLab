"""
@author: Zhou QiYang
@file: run.py
@time: 2025/12/27 15:44
@description: 获取单个实验功能的单元测试
"""

import swanlab


def example_code():
    api = swanlab.OpenApi()
    exp = api.run(path="username/project/expid")  # 可通过api.runs()获取expid
    print(exp.__dict__)
