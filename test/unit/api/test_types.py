#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/14 22:25
@File: test_types.py
@IDE: pycharm
@Description:
    测试开放API模型定义
"""

from swanlab.api.types import Experiment


def test_dict_style_access():
    """
    测试字典风格访问对象字段
    """
    exp = Experiment.model_validate({
        "cuid": "test_cuid",
        "name": "test_experiment",
        "description": "test_description",
        "state": "FINISHED",
        "show": True,
        "createdAt": "2025-05-14T22:25:00Z",
        "finishedAt": "2025-05-14T22:30:00Z",
        "user": {
            "username": "test_user",
            "name": "Test User"
        },
        "profile": {
            "requirements": "test_requirements",
        }
    })
    assert exp["cuid"] == "test_cuid"
    