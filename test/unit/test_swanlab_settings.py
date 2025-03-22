#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-20
@File: test/unit/test_swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab settings module unit tests
"""

import pytest
import swanlab
from swanlab import merge_settings
from swanlab.swanlab_settings import Settings, settings
from pydantic import ValidationError


class TestSwanlabSettingsBasics:
    def test_settings_frozen(self):
        """测试Settings对象实例化后是不可变的"""
        settings = Settings()
        with pytest.raises(ValidationError):
            settings.hardware_monitor = False

    def test_settings_initialization(self):
        """测试Settings的初始化和验证"""
        # 测试基本初始化
        settings = Settings()
        assert isinstance(settings.hardware_monitor, bool)

        # 测试自定义参数
        settings = Settings(hardware_monitor=False)
        assert settings.hardware_monitor is False

        # 测试不允许未定义的字段
        with pytest.raises(ValueError):
            Settings(hardware_monitor=True, unknown_field="value")

    def test_settings_class(self):
        """测试Settings类的基本功能"""
        # 测试默认值
        settings = Settings()
        assert settings.hardware_monitor is True

        # 测试自定义值
        settings = Settings(hardware_monitor=False)
        assert settings.hardware_monitor is False

        # 测试不允许未定义的字段
        with pytest.raises(ValueError):
            Settings(unknown_field="value")


class TestSwanlabSettings:
    def test_settings_creation(self):
        """测试Settings对象的创建"""
        # 测试默认值
        default_settings = Settings()
        assert default_settings.hardware_monitor is True

        # 测试自定义值
        custom_settings = Settings(hardware_monitor=False)
        assert custom_settings.hardware_monitor is False

    def test_merge_settings(self):
        """测试合并设置到全局状态"""
        # 合并自定义设置
        settings = Settings(hardware_monitor=False)
        result = merge_settings(settings)

        # 验证设置已被更新
        assert settings.hardware_monitor is False

        # 再次修改设置
        new_settings = Settings(hardware_monitor=True)
        result = merge_settings(new_settings)

        # 验证设置已被更新
        assert settings.hardware_monitor is True

    def test_default_setup(self):
        """测试不提供设置时的默认行为"""
        # 不提供设置执行init
        swanlab.init()

        # 验证使用了默认设置
        assert settings.hardware_monitor is True

    def test_type_validation(self):
        """测试类型验证"""
        # 传入非Settings对象到merge_settings
        with pytest.raises(TypeError):
            merge_settings({"hardware_monitor": False})

        # 传入非Settings对象到init
        with pytest.raises(TypeError):
            swanlab.init(settings={"hardware_monitor": False})
