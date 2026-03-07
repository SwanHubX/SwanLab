"""
@author: cunyue
@file: test_settings_compatible.py
@time: 2026/3/5 20:59
@description: 测试 SwanLab 一些环境配置的向下兼容提取
"""

from pathlib import Path

from swanlab.sdk.internal.settings import Settings


def test_e2e_legacy_env_compatibility(monkeypatch):
    """测试端到端：旧版扁平环境变量的向下兼容提取"""
    # 模拟旧版环境变量环境
    envs = {
        "SWANLAB_EXP_NAME": "legacy_exp",
        "SWANLAB_DESCRIPTION": "legacy_desc",
        "SWANLAB_TAGS": "tag1, tag2",
        "SWANLAB_WEBHOOK": "https://hook.com",
        "SWANLAB_WEBHOOK_VALUE": "ping",
        "SWANLAB_DASHBOARD_HOST": "0.0.0.0",
        "SWANLAB_DASHBOARD_PORT": "9090",
        "SWANLAB_SAVE_DIR": "/tmp/swanlab",
    }
    for k, v in envs.items():
        monkeypatch.setenv(k, v)

    # 核心：通过根 Settings 实例进行端到端加载
    s = Settings()
    # 验证 root 的兼容映射
    assert s.root == Path("/tmp/swanlab")

    # 验证 experiment 子模块的兼容映射
    assert s.experiment.name == "legacy_exp"
    assert s.experiment.description == "legacy_desc"
    assert s.experiment.tags == ["tag1", "tag2"]

    # 验证 integration 子模块的兼容映射
    assert s.integration.webhook.url == "https://hook.com"
    assert s.integration.webhook.value == "ping"
    assert s.integration.dashboard.host == "0.0.0.0"
    assert s.integration.dashboard.port == 9090


def test_e2e_explicit_input_priority(monkeypatch):
    """测试端到端优先级：代码显式传参 > 环境变量兜底"""
    monkeypatch.setenv("SWANLAB_EXP_NAME", "env_exp")
    monkeypatch.setenv("SWANLAB_WEBHOOK", "env_hook")

    # 模拟在初始化 Settings 时，通过字典/ kwargs 显式传入高优配置
    s = Settings.model_validate(
        {"experiment": {"name": "explicit_exp"}, "integration": {"webhook": {"url": "explicit_hook"}}}
    )
    # 断言：用户的显式输入必须战胜旧版环境变量
    assert s.experiment.name == "explicit_exp"
    assert s.integration.webhook.url == "explicit_hook"
