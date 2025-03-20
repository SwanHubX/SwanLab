#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20
@File: swanlab\swanlab_settings.py
@IDE: vscode
@Description:
    SwanLab全局功能开关，用于管理和控制SwanLab的全局设置
"""

from typing import Dict, Any

from pydantic import BaseModel


class SettingsState:
    """
    单例模式，管理设置状态。
    用于确保全局只有一个设置实例，并提供设置的锁定机制。
    """

    _instance = None
    _initialized = False

    def __new__(cls) -> 'SettingsState':
        if cls._instance is None:
            cls._instance = super(SettingsState, cls).__new__(cls)
            cls._instance._locked = False
            cls._instance._settings: Dict[str, Any] = {}
        return cls._instance

    def lock(self) -> None:
        """锁定设置，防止后续修改"""
        self._locked = True

    @property
    def is_locked(self) -> bool:
        """
        返回设置是否已锁定
        :return: bool 锁定状态
        """
        return self._locked

    def get_settings(self) -> Dict[str, Any]:
        """
        获取当前设置
        :return: Dict 当前设置的副本
        """
        return self._settings.copy()

    def update_settings(self, settings: 'Settings') -> None:
        """
        更新设置，如果已锁定则抛出异常
        :param settings: Settings对象，包含新的设置值
        :raises RuntimeError: 当设置已被锁定时抛出
        """
        if self._locked:
            raise RuntimeError("Cannot modify settings after swanlab.init() has been called")
        self._settings.update(settings.__dict__)


class Settings(BaseModel):
    """
    用户设置类，定义SwanLab支持的所有设置项
    """

    hardware_monitor: bool = True

    def __init__(self, hardware_monitor: bool = True, **kwargs):
        """
        初始化设置
        :param hardware_monitor: 是否启用硬件监控
        :param kwargs: 其他自定义设置项
        """
        self.hardware_monitor = hardware_monitor
        for key, value in kwargs.items():
            setattr(self, key, value)


def merge_settings(settings: Settings) -> Dict[str, Any]:
    """
    合并用户设置到全局设置
    :param settings: Settings对象
    :return: Dict 合并后的设置
    :raises TypeError: 当输入不是Settings对象时抛出
    :raises RuntimeError: 当设置已被锁定时抛出
    """
    if not isinstance(settings, Settings):
        raise TypeError("Expected Settings object")

    state = SettingsState()
    if state.is_locked:
        raise RuntimeError("Cannot modify settings after swanlab.init() has been called")

    state.update_settings(settings)
    return state.get_settings()


def setup(settings: Settings = None) -> Dict[str, Any]:
    """
    初始化单例SettingsState并锁定设置
    :param settings: 可选的Settings对象
    :return: Dict 当前的设置配置
    :raises TypeError: 当settings不是Settings对象时抛出
    """
    state = SettingsState()

    if settings is not None:
        if not isinstance(settings, Settings):
            raise TypeError("Expected Settings object")
        state.update_settings(settings)

    state.lock()

    current_settings = state.get_settings()
    print(f"SwanLab initialized with settings: {current_settings}")

    return current_settings


# 创建默认设置实例
settings = Settings()
