#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20:13
@File: swanlab\swanlab_settings.py
@IDE: pycharm
@Description:
    swanlab全局功能开关
"""

class SettingsState:
    """单例模式，管理设置状态"""
    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(SettingsState, cls).__new__(cls)
            cls._instance._locked = False
            cls._instance._settings = {}
        return cls._instance

    def lock(self):
        """锁定设置，防止后续修改"""
        self._locked = True

    @property
    def is_locked(self):
        """返回设置是否已锁定"""
        return self._locked

    def get_settings(self):
        """获取当前设置"""
        return self._settings.copy()  # 返回副本防止直接修改

    def update_settings(self, settings):
        """更新设置，如果已锁定则抛出异常"""
        if self._locked:
            raise RuntimeError("Cannot modify settings after swanlab.init() has been called")
        self._settings.update(settings.__dict__)


class Settings:
    """用户设置类"""

    def __init__(self, hardware_monitor=True, **kwargs):
        self.hardware_monitor = hardware_monitor
        for key, value in kwargs.items():
            setattr(self, key, value)


def merge_settings(settings):
    """合并用户设置到全局设置"""
    if not isinstance(settings, Settings):
        raise TypeError("Expected Settings object")

    state = SettingsState()
    if state.is_locked:
        raise RuntimeError("Cannot modify settings after swanlab.init() has been called")

    state.update_settings(settings)
    return state.get_settings()


def init(settings=None):
    """初始化SwanLab，并锁定设置"""
    state = SettingsState()

    # 如果传入了新设置，先应用它们
    if settings is not None:
        if not isinstance(settings, Settings):
            raise TypeError("Expected Settings object")
        state.update_settings(settings)

    # 锁定设置，防止后续修改
    state.lock()

    # 使用锁定后的设置进行实际初始化
    current_settings = state.get_settings()
    print(f"SwanLab initialized with settings: {current_settings}")

    # 这里执行真正的初始化逻辑
    # ...

    return current_settings