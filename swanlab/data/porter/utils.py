"""
@author: cunyue
@file: filter.py
@time: 2025/7/20 16:08
@description: 工具函数
"""

from typing import Optional

from swanlab.data.store import RemoteMetric


# 将逻辑写出来的原因是方便测试
def filter_column(key: str, metrics: RemoteMetric) -> bool:
    """
    筛选函数，检查给定的 key 是否在字典中
    :param key: 要检查的 key
    :param metrics: 可选的存储字典，如果为 None，则默认使用 DataPorter 的 run_store.metrics
    :return: 如果 key 在 store 中，则返回 True，否则返回 False
    """
    # 不存在则保留
    if key not in (metrics or {}):
        return True
    return False


# 将逻辑写出来的原因是方便测试
def filter_metric(key: str, step: int, metrics: RemoteMetric) -> bool:
    """
    筛选函数，检查给定的 key 和 step 是否在字典中
    :param key: 要检查的 key
    :param step: 要检查的 step
    :param metrics: 可选的存储字典，如果为 None，则默认使用 DataPorter 的 run_store.metrics
    :return: 如果 key 在 store 中且 step 大于最新 step，则返回 True，否则返回 False
    """
    metrics = metrics or {}
    if key not in metrics:
        return True
    _, _, _, latest_step = metrics[key]
    if step > latest_step:
        return True
    return False


# 将逻辑写出来的原因是方便测试
def filter_epoch(epoch: int, now_epoch: Optional[int]) -> bool:
    """
    筛选函数，检查给定的 epoch 是否大于 now_epoch
    :param epoch: 要检查的 epoch
    :param now_epoch: 当前的 epoch
    :return: 如果 epoch 大于 now_epoch，则返回 True，否则返回 False
    """
    return epoch > (now_epoch or -1)
