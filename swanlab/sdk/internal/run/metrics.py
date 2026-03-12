"""
@author: cunyue
@file: metrics.py
@time: 2026/3/12 14:26
@description: SwanLab 运行时指标管理
"""

from typing import Optional

from swanlab.sdk.internal.context import RunContext


def next_step(ctx: RunContext, user_step: Optional[int] = None) -> int:
    """
    获取下一个全局步数
    """
    metrics = ctx.metrics
    with metrics.lock:
        if user_step is not None:
            metrics.global_step = max(metrics.global_step, user_step)
            return user_step

        metrics.global_step += 1
        return metrics.global_step


def has_metric(ctx: RunContext, key: str) -> bool:
    """
    检查是否存在指定的指标
    """
    with ctx.metrics.lock:
        return ctx.metrics.has_metric(key)
