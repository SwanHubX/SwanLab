"""
@author: cunyue
@file: cloud.py
@time: 2026/3/11 12:15
@description: SwanLab云回调
"""

from swanlab.sdk.utils.callbacker import SwanLabCallback


class CloudCallback(SwanLabCallback):
    @property
    def name(self) -> str:
        return "__swanlab__.cloud"
