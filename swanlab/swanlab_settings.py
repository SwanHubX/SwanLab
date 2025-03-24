#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20
@File: swanlab/swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab全局功能开关，用于管理和控制SwanLab的全局设置
"""

import math
from typing import Any, Optional

try:
    from typing import Annotated  # Python 3.9+
except ImportError:
    from typing_extensions import Annotated  # Python 3.8

from pydantic import BaseModel, ConfigDict, PositiveFloat, PositiveInt, AfterValidator, StrictBool


def is_power_of_two(value: int) -> int:
    """检查一个数是否是2的幂次方"""

    if value <= 0 or not math.log2(value).is_integer():
        raise ValueError(f'{value} is not a power of 2')
    return value


class Settings(BaseModel):
    """
    用户设置类，定义SwanLab支持的所有设置项
    """

    # 使用 Pydantic 的冻结功能代替 @dataclass(frozen=True)
    model_config = ConfigDict(frozen=True)

    # ---------------------------------- 硬件监控部分 ----------------------------------
    # 是否开启硬件监控，如果元数据获取被关闭，则此项无效
    hardware_monitor: Optional[StrictBool] = None
    # 是否开启 CPU 监控，如果硬件监控关闭，则此项无效
    enable_cpu_monitor: Optional[StrictBool] = None
    # 是否开启 GPU 监控，如果硬件监控关闭，则此项无效
    enable_gpu_monitor: Optional[StrictBool] = None
    # 是否开启内存监控，如果硬件监控关闭，则此项无效
    enable_memory_monitor: Optional[StrictBool] = None

    # ---------------------------------- 日志上传部分 ----------------------------------
    # 日志上传间隔
    upload_interval: Optional[PositiveFloat] = None
    # 终端日志上传单行最大字符数
    max_log_line_length: Optional[PositiveInt] = None

    # 主要用来演示 validator 的用法
    memory_block_size: Optional[Annotated[int, AfterValidator(is_power_of_two)]] = None

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


settings = Settings()


def get_settings():
    """获取当前全局设置"""
    global settings
    return settings


def set_settings(new_settings: Settings) -> Settings:
    global settings
    settings = new_settings
    return settings


def reset_settings():
    global settings
    settings = Settings()
