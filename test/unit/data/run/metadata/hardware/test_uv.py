#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@author: cunyue
@DATE: 2025-04-09 14:48:33
@File: test/unit/data/run/metadata/hardware/test_uv.py
@IDE: vscode
@Description:
    测试uv
"""
import subprocess

import pytest
import yaml

from swanlab.data.run.metadata.uv import get_uv


has_uv = subprocess.run(["uv", "--version"], capture_output=True).returncode == 0


@pytest.mark.skipif(not has_uv, reason="uv is not installed")
def test_uv():
    uv_info = get_uv()
    assert isinstance(uv_info, str)
