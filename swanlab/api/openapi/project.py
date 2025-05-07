#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: project.py
@IDE: pycharm
@Description:
    项目相关的开放API
"""

from swanlab.api.http import HTTP
from swanlab.api.openapi.base import ApiBase, ApiHTTP


class ProjectAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)
