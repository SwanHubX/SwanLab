#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20
@File: swanlab/swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab全局功能开关，用于管理和控制SwanLab的全局设置
"""

from typing import Any, Annotated
from pydantic import BaseModel, ConfigDict, PositiveFloat, PositiveInt, AfterValidator


def is_power_of_two(value: int) -> bool:
    """检查一个数是否是2的幂次方"""
    import math

    if value <= 0 or not math.log2(value).is_integer():
        raise ValueError(f'{value} is not a power of 2')
    return value


class Settings(BaseModel):
    """
    用户设置类，定义SwanLab支持的所有设置项

    Attributes
    ----------
    hardware_monitor : bool, optional
        是否启用硬件监控功能。当设置为True时，SwanLab将收集并记录系统硬件信息。
        默认为None，表示使用系统默认设置。

    log_upload : bool, optional
        是否启用日志上传功能。当设置为True时，SwanLab将自动上传日志到云端服务器。
        默认为None，表示使用系统默认设置。

    upload_interval : PositiveFloat, optional
        上传线程的上传间隔，单位为秒。此参数控制SwanLab上传数据到云端的频率。
        必须为正数。默认为None，表示使用系统默认上传间隔。

    enable_cpu_monitor : bool, optional
        是否启用CPU监控。当设置为True时，SwanLab将收集并记录CPU使用率、温度等信息。
        默认为None，表示使用系统默认设置。

    enable_gpu_monitor : bool, optional
        是否启用GPU监控。当设置为True时，SwanLab将收集并记录GPU使用率、显存等信息。
        默认为None，表示使用系统默认设置。

    enable_memory_monitor : bool, optional
        是否启用内存监控。当设置为True时，SwanLab将收集并记录系统内存使用情况。
        默认为None，表示使用系统默认设置。

    max_log_line_length : PositiveInt, optional
        终端实验日志单行最大字符数。超过此长度的日志行将被截断。
        必须为正整数。默认为None，表示使用系统默认值。

    experiment_id : str, optional
        实验ID，用于支持实验断点续训功能。提供此ID可以恢复之前的实验状态。
        默认为None，表示创建新的实验。

    memory_block_size : Annotated[int, AfterValidator(is_power_of_two)], optional
        GPU内存分配的块大小（字节）。必须是2的幂次方（如256, 512, 1024等），
        这有助于优化GPU内存分配和减少内存碎片。如果提供的值不是2的幂次方，
        将引发ValueError。默认为None，表示使用框架默认值。
    """

    # 使用 Pydantic 的冻结功能代替 @dataclass(frozen=True)
    model_config = ConfigDict(frozen=True)

    hardware_monitor: bool | None = None
    log_upload: bool | None = None
    upload_interval: PositiveFloat | None = None
    enable_cpu_monitor: bool | None = None
    enable_gpu_monitor: bool | None = None
    enable_memory_monitor: bool | None = None
    max_log_line_length: PositiveInt | None = None
    experiment_id: str | None = None

    # 主要用来演示 validator 的用法
    memory_block_size: Annotated[int, AfterValidator(is_power_of_two)] | None = None

    def get(self, field_name: str, default: Any = None) -> Any:
        """
        获取设置项的值，当值为None时返回默认值
        :param field_name: 字段名称
        :param default: 默认值，当字段值为None时返回
        :return: 字段值或默认值
        :raises KeyError: 当字段名不存在时抛出
        """
        if field_name not in self.model_fields:
            raise KeyError(f"Field '{field_name}' does not exist")

        value = getattr(self, field_name)
        return default if value is None else value


settings = Settings(
    hardware_monitor=True,
    log_upload=True,
    upload_interval=10,
    enable_cpu_monitor=True,
    enable_gpu_monitor=True,
    enable_memory_monitor=True,
    max_log_line_length=100,
    experiment_id=None,
    memory_block_size=1024,
)
