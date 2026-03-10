"""
@author: cunyue
@file: test_settings_integration.py
@time: 2026/3/10 15:32
@description: 测试 SwanLab 集成配置
"""

from swanlab.sdk.internal.settings import Settings


def test_integration_nested_env_parse(monkeypatch):
    """
    测试通过 Pydantic 标准嵌套环境变量解析 Integration 配置，
    验证 max_split=1 的拦截器是否生效。
    """
    monkeypatch.setenv("SWANLAB_INTEGRATION_WEBHOOK_URL", "https://nested.hook.com")
    monkeypatch.setenv("SWANLAB_INTEGRATION_WEBHOOK_VALUE", "nested_token")
    monkeypatch.setenv("SWANLAB_INTEGRATION_DASHBOARD_HOST", "192.168.1.100")
    monkeypatch.setenv("SWANLAB_INTEGRATION_DASHBOARD_PORT", "9090")

    s = Settings()

    assert s.integration.webhook.url == "https://nested.hook.com"
    assert s.integration.webhook.value == "nested_token"
    assert s.integration.dashboard.host == "192.168.1.100"
    assert s.integration.dashboard.port == 9090


def test_integration_env_priority(monkeypatch):
    """
    测试当同时存在新旧环境变量时，新版嵌套变量（通过 dict 重组）应当覆盖 default_factory 的旧变量
    """
    monkeypatch.setenv("SWANLAB_WEBHOOK", "https://legacy.hook.com")  # 应该被覆盖
    monkeypatch.setenv("SWANLAB_INTEGRATION_WEBHOOK_URL", "https://new.hook.com")  # 优先级更高

    s = Settings()

    # 验证新版变量生效
    assert s.integration.webhook.url == "https://new.hook.com"
