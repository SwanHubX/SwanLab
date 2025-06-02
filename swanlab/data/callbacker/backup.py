"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

from swankit.callback import SwanKitCallback


class BackupCallback(SwanKitCallback):
    def __str__(self) -> str:
        return "SwanLabBackupCallback"
