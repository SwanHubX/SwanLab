#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/5 16:30
@File: error.py
@IDE: pycharm
@Description:
    错误模型
"""


class OperateErrorInfo:
    """
    操作错误信息，当操作员操作回调发生错误时，会生成这个对象，传给相应的回调函数
    """

    def __init__(self, desc: str):
        self.desc = desc
        """
        错误描述
        """

    def __str__(self):
        return f"SwanLabError: {self.desc}"
