"""
@author: cunyue
@file: swanlab.py
@time: 2026/4/24 12:52
@description: SwanLab 信息采集模块
"""

from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.internal.pkg.helper import get_swanlab_version
from swanlab.sdk.internal.probe_python.context import ProbeContext
from swanlab.sdk.typings.probe_python import SwanLabSnapshot


@safe.decorator(level="debug", message="Failed to get swanlab environment")
def get(ctx: ProbeContext) -> SwanLabSnapshot:
    """获取 SwanLab 信息快照"""
    return SwanLabSnapshot(version=get_swanlab_version(), run_dir=str(ctx.config.run_dir))
