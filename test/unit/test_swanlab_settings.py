#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-20
@File: test/unit/test_swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab settings module unit tests
"""
import platform

import pytest
from pydantic import ValidationError

import swanlab
from swanlab.swanlab_settings import get_settings, reset_settings
from tutils import TEMP_PATH


class TestSwanlabSettingsBasics:
    def test_settings_frozen(self):
        """测试Settings对象实例化后是不可变的"""
        settings = swanlab.Settings()
        with pytest.raises(ValidationError):
            settings.hardware_monitor = False  # noqa

    def test_settings_initialization(self):
        """测试Settings的初始化和验证"""
        # 测试基本初始化
        settings = swanlab.Settings()
        assert settings.hardware_monitor is True

        # 测试自定义参数
        settings = swanlab.Settings(hardware_monitor=False)
        assert settings.hardware_monitor is False


def test_log_proxy_type():
    """测试日志代理类型"""
    # 测试默认值
    settings = swanlab.Settings()
    assert settings.log_proxy_type == "all"

    # 测试自定义值
    settings = swanlab.Settings(log_proxy_type="stdout")
    assert settings.log_proxy_type == "stdout"

    # 测试无效值
    with pytest.raises(ValidationError):
        swanlab.Settings(log_proxy_type="invalid")


class TestSwanlabSettings:
    def teardown_method(self):
        # 每个测试方法后执行
        try:
            swanlab.finish()
        except Exception:
            pass
        reset_settings()

    def test_merge_settings(self):
        """测试合并设置到全局状态"""
        # 合并自定义设置
        settings1 = swanlab.Settings(hardware_monitor=False)
        swanlab.merge_settings(settings1)

        settings = get_settings()
        assert settings.hardware_monitor is False

        # 再次修改设置
        settings2 = swanlab.Settings(hardware_monitor=True)
        swanlab.merge_settings(settings2)

        settings = get_settings()
        assert settings.hardware_monitor is True

    def test_default_setup(self):
        """测试不提供设置时的默认行为"""
        # 不提供设置执行init
        swanlab.init(mode="disabled")

        # 验证使用了默认设置
        settings = get_settings()
        assert settings.hardware_monitor is True

    def test_change_settings(self):
        """测试修改设置"""
        settings = swanlab.Settings(hardware_monitor=False)
        swanlab.merge_settings(settings)

        settings = get_settings()
        assert settings.hardware_monitor is False

        new_settings = swanlab.Settings(hardware_monitor=True)
        swanlab.merge_settings(new_settings)

        settings = get_settings()
        assert settings.hardware_monitor is True

    def test_change_settings_with_init(self):
        """测试在init时修改设置"""
        settings = swanlab.Settings(hardware_monitor=False)
        swanlab.merge_settings(settings)
        settings = get_settings()  # 此时硬件监控被关闭
        assert settings.hardware_monitor is False

        new_settings = swanlab.Settings(hardware_monitor=True)
        swanlab.init(
            settings=new_settings,
            mode="disabled",
        )  # 此时硬件监控被开启
        settings = get_settings()
        assert settings.hardware_monitor is True

    def test_type_validation(self):
        """测试类型验证"""
        # 传入非Settings对象到merge_settings
        with pytest.raises(TypeError):
            swanlab.merge_settings({"hardware_monitor": False})

        # 传入非Settings对象到init
        with pytest.raises(TypeError):
            swanlab.init(settings={"hardware_monitor": False})

        # 传入基本类型到merge_settings
        with pytest.raises(TypeError):
            swanlab.merge_settings(42)

        with pytest.raises(TypeError):
            swanlab.merge_settings("settings")

        with pytest.raises(TypeError):
            swanlab.merge_settings(True)

        # 传入其他非Settings类的对象
        class FakeSettings:
            hardware_monitor = False

        with pytest.raises(TypeError):
            swanlab.merge_settings(FakeSettings())

        # 测试嵌套对象类型验证
        with pytest.raises(TypeError):
            swanlab.init(settings={"hardware_monitor": False, "nested": {"invalid": True}})

    def test_max_log_length(self):
        """
        测试日志长度范围
        """
        swanlab.Settings(max_log_length=500)
        swanlab.Settings(max_log_length=4096)
        with pytest.raises(ValidationError):
            swanlab.Settings(max_log_length=499)

        with pytest.raises(ValidationError):
            swanlab.Settings(max_log_length=4097)

    def test_disk_io_default_dir(self):
        """
        测试默认磁盘IO监控路径
        """
        settings = swanlab.Settings()
        assert settings.disk_io_dir == "C:\\" if platform.system() == 'Windows' else "/"

    def test_hardware_interval(self):
        settings = swanlab.Settings()
        assert settings.hardware_interval is None
        # 最小值为5
        with pytest.raises(ValidationError):
            swanlab.Settings(hardware_interval=1)
        # 不允许为小数
        with pytest.raises(ValidationError):
            swanlab.Settings(hardware_interval=5.5)


def test_filter_changed_fields():
    # 默认
    settings = swanlab.Settings()
    changed = settings.filter_changed_fields()
    assert changed == {}
    # 修改
    settings = swanlab.Settings(hardware_monitor=False)
    changed = settings.filter_changed_fields()
    assert changed == {"hardware_monitor": False}
    # 修改，但是没有变化
    settings = swanlab.Settings(hardware_monitor=True)
    changed = settings.filter_changed_fields()
    assert changed == {}

    # 修改一些特殊不可被JSON序列化的字段
    settings = swanlab.Settings(disk_io_dir=TEMP_PATH)
    changed = settings.filter_changed_fields()
    assert changed == {"disk_io_dir": TEMP_PATH}


@pytest.mark.parametrize("mode", ["disabled", 'offline', 'cloud', 'local'])
def test_app_settings_ok(mode):
    settings = swanlab.swanlab_settings.FolderSettings(mode=mode)
    assert settings.mode == mode


def test_app_settings_invalid():
    """错误的 mode 会被重置为 cloud，其他参数会被忽略"""
    settings = swanlab.swanlab_settings.FolderSettings(mode="123", x="y")
    # 默认为 cloud
    assert settings.mode == "cloud"


def test_read_folder_settings_ok(tmp_path):
    """读取存在且格式正确的配置文件"""
    file_path = tmp_path / "settings.ini"
    with open(file_path, "w") as f:
        f.write("[default]\nmode=offline\n")
    settings = swanlab.swanlab_settings.read_folder_settings(folder_path=tmp_path, filename='settings.ini')
    assert settings.mode == "offline"


def test_read_folder_settings_invalid(tmp_path):
    """读取存在但是格式错误的配置文件"""
    file_path = tmp_path / "settings.ini"
    with open(file_path, "w") as f:
        f.write("invalid content")
    settings = swanlab.swanlab_settings.read_folder_settings(folder_path=tmp_path, filename='settings.ini')
    assert settings.mode == "cloud"


def test_read_folder_settings_not_exist(tmp_path):
    """
    - 对于不存在的文件不会报错
    - 对于存在但是格式错误的文件也不会报错
    """
    settings = swanlab.swanlab_settings.read_folder_settings(folder_path=tmp_path, filename="not_exist.ini")
    assert settings.mode == "cloud"
    open(tmp_path / "not_exist.ini", "w").close()
    settings = swanlab.swanlab_settings.read_folder_settings(folder_path=tmp_path, filename="not_exist.ini")
    assert settings.mode == "cloud"


def test_write_folder_settings_ok(tmp_path):
    """测试写入配置文件"""
    file_path = tmp_path / "settings.ini"
    swanlab.swanlab_settings.write_folder_settings(tmp_path, {'mode': 'local'}, filename="settings.ini")
    assert file_path.exists()
    with open(file_path, "r") as f:
        content = f.read()
    assert "[default]" in content
    assert "mode = local" in content
    # 测试读取
    settings = swanlab.swanlab_settings.read_folder_settings(folder_path=tmp_path, filename="settings.ini")
    assert settings.mode == "local"
