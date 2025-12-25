"""
@author: Zhou Qiyang
@file: runs.py
@time: 2025/12/24 16:58
@description: 
"""

import random
from unittest.mock import patch

import swanlab

# 测试实验数据
test_runs = [
    [
        [
            {
                "job": None,
                "cuid": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=20)),
                "name": f"run-test-{n}",
                "show": True,
                "pin": False,
                "state": ["FINISHED", "RUNNING", "CRASHED", "ABORTED"][n % 4],
                "colors": ["#b15fbb", "#b15fbb"],
                "cluster": "GG",
                "cloneType": "COPY",
                "rootProId": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=10)),
                "rootExpId": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=10)),
                "createdAt": "2025-01-01T01:00:00Z",
                "updatedAt": "2025-01-01T01:00:00Z",
                "finishedAt": "2025-01-01T03:00:00Z",
                "pinnedAt": None,
                "description": "fake exp",
                "labels": [{"name": "label1"}, {"name": "label2"}],
                "user": {"username": "username", "name": "name"},
                "profile": {"config": {"lr": 0.00001}, "scalar": {"acc": 0.92}},
                "resumes": [],
            }
            for n in range(50)
        ]
    ]
]

# 分组时的实验数据
grouped_test_runs = [[{"GG": {"user": test_runs[0][0]}}]]


def test_get_runs():
    """测试获取所有实验的基本功能"""
    with patch("swanlab.core_python.client.Client.post", side_effect=test_runs):
        api = swanlab.OpenApi()
        exps = api.runs(path="user/test-project")

        for run, index in zip(exps, range(50)):
            # 测试基本属性
            assert run.name == f"run-test-{index}"
            assert run.state == ["FINISHED", "RUNNING", "CRASHED", "ABORTED"][index % 4]
            assert run.description == "fake exp"
            assert run.group == "GG"
            assert run.history_line_count == 50
            assert run.job is None
            assert run.created_at == "2025-01-01T01:00:00Z"

            # 测试用户信息
            assert run.user.username == "username"
            assert run.user.name == "name"

            # 测试标签
            assert len(run.labels) == 2
            assert run.labels[0].name == "label1"
            assert run.labels[1].name == "label2"

            # 测试配置和摘要
            assert run.config == {"lr": 0.00001}
            assert run.summary == {"acc": 0.92}

            # 测试其他属性
            assert run.metric_keys == ["acc"]
            assert run.history_line_count == 50


def test_get_grouped_runs():
    """测试分组情况下，返回格式字典的实验数据"""
    with patch("swanlab.core_python.client.Client.post", side_effect=grouped_test_runs):
        api = swanlab.OpenApi()
        exps = api.runs(path="user/test-project")

        for run, index in zip(exps, range(50)):
            assert run.name == f"run-test-{index}"
            assert run.history_line_count == 50
