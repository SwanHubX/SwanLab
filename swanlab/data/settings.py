#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:45:15
@File: swanlab/database/settings.py
@IDE: vscode
@Description:
    数据收集部分配置，此为运行时生成的配置，
"""
import os
from ..env import root


class SwanDataSettings:
    def __init__(self, exp_name: str) -> None:
        """实验名称

        Parameters
        ----------
        exp_name : str
            _description_
        """
        self.exp_name: str = exp_name
