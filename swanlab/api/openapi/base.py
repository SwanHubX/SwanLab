#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:36
@File: base.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI API基类
"""

from swanlab.api.http import HTTP


class ApiBase:
    def __init__(self, http: HTTP):
        self.http: HTTP = http
