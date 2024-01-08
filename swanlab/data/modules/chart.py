#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-29 12:25:51
@File: swanlab/database/modules/chart.py
@IDE: vscode
@Description:
    SwanLab图表选择器，方便用户选择受支持的图表类型
"""


class Chart:
    default = "line", [float, int]
    """默认图表，为折线图，允许接收`float, int`"""
    line = "line", [float, int]
    """折线图，允许接收`float, int`"""
    image = "image", [str]
    """图片，允许接收`str`"""
