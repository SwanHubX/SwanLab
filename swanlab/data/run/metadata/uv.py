#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@author: KashiwaByte
@DATE: 2025-04-09 13:57:37
@File: swanlab/data/run/metadata/uv.py
@IDE: vscode
@Description:
    收集uv环境信息
"""
import subprocess
from typing import Optional

import yaml


def get_uv() -> Optional[str]:
    """
    获取uv信息，如果不存在则返回None
    Returns: str
    """
    result = subprocess.run(["uv", "pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=15)
    # 检查命令是否成功运行
    if result.returncode == 0:
        return result.stdout
    else:
        return None
