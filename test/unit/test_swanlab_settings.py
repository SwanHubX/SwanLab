#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-20
@File: test\test_swanlab_settings.py
@IDE: vscode
@Description:
    SwanLab settings module unit tests
"""

import pytest
from swanlab.swanlab_settings import SettingsState, Settings, merge_settings, setup


def test_settings_state_singleton():
    """测试SettingsState的单例模式"""
    # 创建两个实例，应该是同一个对象
    state1 = SettingsState()
    state2 = SettingsState()
    assert state1 is state2

    # 确保设置是共享的
    state1._settings['test'] = 'value'
    assert state2._settings['test'] == 'value'


def test_settings_state_lock():
    """测试SettingsState的锁定机制"""
    state = SettingsState()
    state._settings.clear()  # 清除之前测试可能留下的设置

    # 测试锁定前可以更新设置
    settings = Settings(hardware_monitor=True)
    state.update_settings(settings)
    assert state._settings['hardware_monitor'] is True

    # 锁定后不能更新设置
    state.lock()
    assert state.is_locked is True

    with pytest.raises(RuntimeError):
        state.update_settings(settings)


def test_settings_class():
    """测试Settings类的基本功能"""
    # 测试默认值
    settings = Settings()
    assert settings.hardware_monitor is True

    # 测试自定义值
    settings = Settings(hardware_monitor=False)
    assert settings.hardware_monitor is False

    # 测试动态添加设置
    settings = Settings(hardware_monitor=True, custom_setting="test")
    assert settings.hardware_monitor is True
    assert settings.custom_setting == "test"


def test_merge_settings():
    """测试merge_settings函数"""
    # 重置SettingsState
    state = SettingsState()
    state._settings.clear()
    state._locked = False

    # 测试合并有效的设置
    settings = Settings(hardware_monitor=False, custom_setting="test")
    result = merge_settings(settings)
    assert result['hardware_monitor'] is False
    assert result['custom_setting'] == "test"

    # 测试无效输入
    with pytest.raises(TypeError):
        merge_settings("invalid")

    # 测试锁定后的合并
    state.lock()
    with pytest.raises(RuntimeError):
        merge_settings(settings)


def test_setup():
    """测试setup函数"""
    # 重置SettingsState
    state = SettingsState()
    state._settings.clear()
    state._locked = False

    # 测试默认设置
    result = setup()
    assert isinstance(result, dict)
    assert 'hardware_monitor' in result
    assert result['hardware_monitor'] is True

    # 重置状态
    state._settings.clear()
    state._locked = False

    # 测试自定义设置
    custom_settings = Settings(hardware_monitor=False, custom_setting="test")
    result = setup(custom_settings)
    assert result['hardware_monitor'] is False
    assert result['custom_setting'] == "test"

    # 测试无效输入
    with pytest.raises(TypeError):
        setup("invalid")

    # 测试重复setup
    with pytest.raises(RuntimeError):
        setup(custom_settings)


def test_settings_copy():
    """测试设置的复制机制"""
    state = SettingsState()
    state._settings.clear()
    state._locked = False

    # 设置初始值
    settings = Settings(hardware_monitor=True)
    state.update_settings(settings)

    # 获取设置的副本
    settings_copy = state.get_settings()

    # 修改副本不应影响原始设置
    settings_copy['hardware_monitor'] = False
    assert state._settings['hardware_monitor'] is True


def test_settings_initialization():
    """测试Settings的初始化和验证"""
    # 测试基本初始化
    settings = Settings()
    assert isinstance(settings.hardware_monitor, bool)

    # 测试自定义参数
    custom_settings = {'hardware_monitor': False, 'log_level': 'DEBUG', 'max_retries': 3}
    settings = Settings(**custom_settings)
    assert settings.hardware_monitor is False
    assert settings.log_level == 'DEBUG'
    assert settings.max_retries == 3
