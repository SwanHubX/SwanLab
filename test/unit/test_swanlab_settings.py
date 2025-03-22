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

# from swanlab import swanlab_settings
from pydantic import ValidationError


class TestSwanlabSettingsBasics:
    def test_settings_frozen(self):
        """测试Settings对象实例化后是不可变的"""
        settings = swanlab.Settings()
        with pytest.raises(ValidationError):
            settings.hardware_monitor = False

    def test_settings_initialization(self):
        """测试Settings的初始化和验证"""
        # 测试基本初始化
        settings = swanlab.Settings()
        assert settings.hardware_monitor is None

        # 测试自定义参数
        settings = swanlab.Settings(hardware_monitor=False)
        assert settings.hardware_monitor is False

    def test_settings_get(self):
        """测试Settings对象的get方法"""
        settings = swanlab.Settings()
        assert settings.get("hardware_monitor") is None
        assert settings.get("hardware_monitor", 1) == 1


class TestSwanlabSettings:
    def teardown_method(self):
        # 每个测试方法后执行
        try:
            swanlab.finish()
        except Exception:
            pass

    def test_merge_settings(self):
        """测试合并设置到全局状态"""
        # 合并自定义设置
        settings1 = swanlab.Settings(hardware_monitor=False)
        swanlab.merge_settings(settings1)

        settings = swanlab.get_settings()
        assert settings.hardware_monitor is False

        # 再次修改设置
        settings2 = swanlab.Settings(hardware_monitor=True)
        swanlab.merge_settings(settings2)

        settings = swanlab.get_settings()
        assert settings.hardware_monitor is True

    def test_default_setup(self):
        """测试不提供设置时的默认行为"""
        # 不提供设置执行init
        swanlab.init(mode="disabled")

        # 验证使用了默认设置
        settings = swanlab.get_settings()
        assert settings.hardware_monitor is True

    def test_type_validation(self):
        """测试类型验证"""
        # 传入非Settings对象到merge_settings
        with pytest.raises(TypeError):
            swanlab.merge_settings({"hardware_monitor": False})

        # 传入非Settings对象到init
        with pytest.raises(TypeError):
            swanlab.init(settings={"hardware_monitor": False})
