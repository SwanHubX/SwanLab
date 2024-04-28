#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-26 17:08:07
@File: swanlab\integration\zhipuai\zhipuai.py
@IDE: vscode
@Description:
    ZhipuAI autolog 集成
"""
import os
import sys
import zhipuai

current_dir = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.join(current_dir, "..")
sys.path.append(relative_path)

from integration_utils.autologging import AutologAPI

from swanlab.integration.zhipuai.resolver import ZhipuAIClientResponseResolver
import pkg_resources


def get_library_version(library_name):
    try:
        version = pkg_resources.get_distribution(library_name).version
        return version
    except pkg_resources.DistributionNotFound:
        return "Library not found"


# 调用函数并传入库名称
library_name = "zhipuai"
version = get_library_version(library_name)


autolog = AutologAPI(
    name="ZhipuAI",
    symbols=(
        "chat.completions.create",
        "chat.asyncCompletions.create",
    ),
    resolver=ZhipuAIClientResponseResolver(),
    client=zhipuai.ZhipuAI(),
    lib_version=version,
)
