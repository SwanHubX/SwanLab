#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-23 11:38:56
@File: swanlab\intergration\openai\openai.py
@IDE: vscode
@Description:
    OpenAI 集成 autolog
"""
import os
import sys
import openai

current_dir = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.join(current_dir, "..")
sys.path.append(relative_path)

from intergration_utils.autologging import AutologAPI

from intergration.openai.resolver import OpenAIRequestResponseResolver
import pkg_resources


def get_library_version(library_name):
    try:
        version = pkg_resources.get_distribution(library_name).version
        return version
    except pkg_resources.DistributionNotFound:
        return "Library not found"


# 调用函数并传入库名称
library_name = "openai"
version = get_library_version(library_name)

if version == "0.28.0":
    autolog = AutologAPI(
        name="OpenAI",
        symbols=(
            # "OpenAI().Chat.Completion.create",
            "ChatCompletion.create",
            "Edit.create",
            "Completion.create",
            "Edit.acreate",
            "Completion.acreate",
            "ChatCompletion.acreate",
        ),
        resolver=OpenAIRequestResponseResolver(),
    )
else:
    # 1.0版本以后的OpenAI
    autolog = AutologAPI(
        name="OpenAI",
        symbols=(
            # "OpenAI().Chat.Completion.create",
            "chat.completions.create",
            # "Edit.create",
            # "Completion.create",
            # "Edit.acreate",
            # "Completion.acreate",
            # "ChatCompletion.acreate",
        ),
        resolver=OpenAIRequestResponseResolver(),
        client=openai.OpenAI(),
    )
