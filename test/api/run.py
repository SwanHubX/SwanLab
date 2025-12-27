"""
@author: Zhou QiYang
@file: run.py
@time: 2025/12/27 15:44
@description: 获取单个实验功能的单元测试
"""

import swanlab


def test_get_single_exp():
    api = swanlab.OpenApi()
    exp = api.run(path="bainiantest/resume_test/f59muc3gkbwkjns37b9aw")
    print(exp.__dict__)
