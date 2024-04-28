#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-23 11:38:56
@File: swanlab\intergration\openai\openai.py
@IDE: vscode
@Description:
    OpenAI 集成 autolog
"""
import openai
from swanlab.integration.integration_utils.autologging import AutologAPI
from swanlab.integration.openai.resolver import OpenAIRequestResponseResolver, OpenAIClientResponseResolver
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

# 判断openai版本
if version[0] == "0":
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
        lib_version=version,
    )
else:
    # 1.0版本以后的OpenAI
    autolog = AutologAPI(
        name="OpenAI",
        symbols=(
            "chat.completions.create",
            "completions.create",
            # "Asynccompletions.create",
            # "chat.completions.acreate",
            # "Edit.acreate",
            # "Edit.create",
        ),
        resolver=OpenAIClientResponseResolver(),
        client=openai.OpenAI(),
        lib_version=version,
    )
