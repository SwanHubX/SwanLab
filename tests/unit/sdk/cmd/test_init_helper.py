"""
@author: cunyue
@file: test_init_helper.py
@time: 2026/3/10 20:52
@description: 测试 init 中的辅助函数
"""

import json
from unittest.mock import MagicMock

import pytest

# 请替换为你的实际导入路径
from swanlab.sdk.cmd.init import load_config, prompt_init_mode


# ==========================================
# Fixtures
# ==========================================
@pytest.fixture
def mock_settings():
    """构造一个轻量级的 Settings Mock 对象"""
    settings = MagicMock()
    settings.run.config = None
    settings.mode = "cloud"
    settings.interactive = True
    settings.web_host = "https://swanlab.cn"
    return settings


# ==========================================
# Tests for `load_config`
# ==========================================
def test_load_config_basic(mock_settings):
    """测试基础类型的加载 (None 和 dict)"""
    assert load_config(mock_settings, None) == {}
    assert load_config(mock_settings, {"key": "value"}) == {"key": "value"}


def test_load_config_files(mock_settings, tmp_path):
    """测试从 JSON 和 YAML 文件中正确加载配置"""
    # 测试 JSON
    json_file = tmp_path / "config.json"
    json_file.write_text(json.dumps({"a": 1}))
    assert load_config(mock_settings, str(json_file)) == {"a": 1}

    # 测试 YAML
    yaml_file = tmp_path / "config.yaml"
    yaml_file.write_text("b: 2\nc: 3")
    assert load_config(mock_settings, yaml_file) == {"b": 2, "c": 3}


def test_load_config_exceptions(mock_settings):
    """测试异常分支：文件不存在和类型错误"""
    with pytest.raises(FileNotFoundError):
        load_config(mock_settings, "not_exist_file.json")

    with pytest.raises(ValueError, match="Invalid config type"):
        load_config(mock_settings, 123)  # type: ignore


# ==========================================
# Tests for `prompt_init_mode`
# ==========================================
def test_prompt_fast_exits(mock_settings, monkeypatch):
    """测试快速跳过的条件：非 cloud 模式、非交互模式、已登录"""
    # 1. 模拟已登录 (client.exists 返回 True)
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: True)
    assert prompt_init_mode(mock_settings) == ("cloud", True)

    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)

    # 2. 模拟非交互模式
    mock_settings.interactive = False
    assert prompt_init_mode(mock_settings) == ("cloud", False)

    # 3. 模拟非云端模式
    mock_settings.interactive = True
    mock_settings.mode = "local"
    assert prompt_init_mode(mock_settings) == ("local", False)


def test_prompt_auto_login(mock_settings, monkeypatch):
    """测试本地已存在 apikey 时的自动登录逻辑"""
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)
    monkeypatch.setattr("swanlab.sdk.cmd.init.apikey.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.init.apikey.get", lambda: "fake-key")

    mock_login = MagicMock()
    monkeypatch.setattr("swanlab.sdk.cmd.init.raw_login", mock_login)

    assert prompt_init_mode(mock_settings) == ("cloud", True)
    mock_login.assert_called_once_with(api_key="fake-key")


@pytest.mark.parametrize(
    "inputs, mock_login_success, expected_mode, expected_success",
    [
        # 直接选 3 -> 切换离线模式
        (["3"], False, "offline", False),
        # 直接选 1 -> 触发交互登录并返回成功
        (["1"], True, "cloud", True),
        # 直接选 2 -> 触发交互登录并返回失败
        (["2"], False, "cloud", False),
        # 乱输一通后选 3 -> 循环容错测试
        (["invalid", "wrong", "3"], False, "offline", False),
    ],
)
def test_prompt_interactive_choices(
    mock_settings, monkeypatch, inputs, mock_login_success, expected_mode, expected_success
):
    """使用参数化和猴子补丁极致压缩终端交互的测试代码"""
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)
    monkeypatch.setattr("swanlab.sdk.cmd.init.apikey.exists", lambda: False)

    # 模拟 interactive_login 返回值
    monkeypatch.setattr("swanlab.sdk.cmd.init.interactive_login", lambda save: mock_login_success)

    # 模拟 input，通过迭代器按顺序弹出输入值
    input_iterator = iter(inputs)
    monkeypatch.setattr("builtins.input", lambda prompt: next(input_iterator))

    assert prompt_init_mode(mock_settings) == (expected_mode, expected_success)
