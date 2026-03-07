"""
@author: cunyue
@file: run.py
@time: 2026/3/6 13:07
@description: SwanLab SDK 实验运行类型定义
"""

from typing import Literal

ResumeType = Literal["must", "allow", "never"]
"""
Run 恢复策略
"""

ModeType = Literal["disabled", "cloud", "local", "offline"]
"""
Run 运行策略
"""
