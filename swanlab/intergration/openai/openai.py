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

current_dir = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.join(current_dir, "..")
sys.path.append(relative_path)

from intergration_utils.autologging import AutologAPI

from intergration.openai.resolver import OpenAIRequestResponseResolver


autolog = AutologAPI(
    name="OpenAI",
    symbols=(
        "Edit.create",
        "Completion.create",
        "ChatCompletion.create",
        "Edit.acreate",
        "Completion.acreate",
        "ChatCompletion.acreate",
    ),
    resolver=OpenAIRequestResponseResolver(),
)
