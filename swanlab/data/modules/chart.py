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
    # 默认图表
    default = "line", [float, int]
    # 折线图
    line = "line", [float, int]
    # 图片类型，list代表一步多图
    image = "image", [list, str]
    # 音频类型，list代表一步多音频
    audio = "audio", [list, str]
    # 文本类型，list代表一步多文本
    text = "text", [list, str]
    # 视频类型 list代表一步多视频
    video = "video", [list, str]
    # 3D点云类型，list代表一步多3D点云
    object3d = "object3d", [list, str]
