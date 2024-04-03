#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:28
@File: config.py.py
@IDE: pycharm
@Description:
    存储一些快捷配置
"""
import os
import json
import nanoid

__test_path = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__)
        )
    ),
    "test"
)

TEMP_PATH = os.path.join(__test_path, "temp")

CONFIG: dict = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.json")))

__all__ = ["TEMP_PATH", "CONFIG", "nanoid"]
