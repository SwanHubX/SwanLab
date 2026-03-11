"""
@author: cunyue
@file: local.py
@time: 2026/3/11 12:29
@description: SwanLab SDK 本地回调器
"""

from swanlab.sdk.utils.callbacker import SwanLabCallback


class LocalCallback(SwanLabCallback):
    @property
    def name(self) -> str:
        return "__swanlab__.local"
