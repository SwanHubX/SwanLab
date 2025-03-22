#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20
@File: swanlab/swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab全局功能开关，用于管理和控制SwanLab的全局设置
"""

from typing import Dict, Any
from pydantic import BaseModel, ConfigDict


class Settings(BaseModel):
    """
    用户设置类，定义SwanLab支持的所有设置项
    """

    # 使用 Pydantic 的冻结功能代替 @dataclass(frozen=True)
    model_config = ConfigDict(frozen=True)

    hardware_monitor: bool | None = None
    log_upload: bool | None = None

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
)
