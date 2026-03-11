"""
@author: cunyue
@file: offline.py
@time: 2026/3/11 12:29
@description: SwanLab SDK 离线回调器
"""

from swanlab.sdk.utils.callbacker import SwanLabCallback


class OfflineCallback(SwanLabCallback):
    @property
    def name(self) -> str:
        return "__swanlab__.offline"
