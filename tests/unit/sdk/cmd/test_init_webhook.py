"""
@author: cunyue
@file: test_init_webhook.py
@time: 2026/3/11 15:52
@description: 测试 webhook 回调
"""

import json
from unittest.mock import MagicMock

import pytest
import responses

from swanlab.sdk.cmd.init import send_webhook


@pytest.fixture
def mock_ctx():
    """创建一个模拟的 RunContext 对象"""
    ctx = MagicMock()
    # 设置基础配置
    ctx.run_dir = "/path/to/run_dir"

    # 设置 settings 相关属性
    settings = ctx.config.settings
    settings.mode = "local"
    settings.integration.webhook.url = "https://example.com/webhook"
    settings.integration.webhook.value = "test_webhook_value"
    settings.integration.webhook.timeout = 5

    return ctx


@pytest.fixture
def setup_mocks(monkeypatch):
    """使用 monkeypatch 模拟外部依赖"""
    # 模拟 get_swanlab_version
    # 请确保 "your_module.get_swanlab_version" 是正确的导入路径
    monkeypatch.setattr("swanlab.sdk.cmd.init.get_swanlab_version", lambda: "1.0.0")

    # 模拟 console.debug
    mock_console = MagicMock()
    monkeypatch.setattr("swanlab.sdk.cmd.init.console", mock_console)

    return mock_console


def test_send_webhook_disabled(mock_ctx, setup_mocks):
    """测试 mode 为 disabled 的情况，断言提早退出并未发送请求"""
    # 动态修改为 disabled
    mock_ctx.config.settings.mode = "disabled"

    # from your_module import send_webhook
    send_webhook(mock_ctx)

    # 验证 debug 信息被正确打印
    setup_mocks.debug.assert_called_once_with("Skipping webhook because mode is disabled.")


def test_send_webhook_no_url(mock_ctx, setup_mocks):
    """测试没有设置 webhook URL 的情况"""
    mock_ctx.config.settings.integration.webhook.url = None

    send_webhook(mock_ctx)

    setup_mocks.debug.assert_called_once_with("Skipping webhook because SWANLAB_WEBHOOK is not set.")


@responses.activate
def test_send_webhook_cloud_mode(mock_ctx, setup_mocks):
    """测试 cloud 模式下的正常请求，验证云端 exp_url 拼接和 JSON 载荷"""
    settings = mock_ctx.config.settings
    settings.mode = "cloud"
    settings.web_host = "https://swanlab.cn"
    settings.project.workspace = "test_user"
    settings.project.name = "test_project"
    settings.run.id = "run_12345"

    webhook_url = settings.integration.webhook.url

    # 注册 responses 拦截
    responses.add(responses.POST, webhook_url, status=200)

    send_webhook(mock_ctx)

    # 验证是否发送了请求
    assert len(responses.calls) == 1
    request = responses.calls[0].request
    assert request.url == webhook_url

    # 验证请求的 JSON 体是否符合预期
    payload = json.loads(request.body)  # type: ignore
    assert payload == {
        "value": "test_webhook_value",
        "swanlab": {
            "version": "1.0.0",
            "mode": "cloud",
            "run_dir": "/path/to/run_dir",
            "exp_url": "https://swanlab.cn/@test_user/test_project/runs/run_12345",
        },
    }


@responses.activate
def test_send_webhook_local_mode(mock_ctx, setup_mocks):
    """测试 local 模式下的正常请求，验证 exp_url 为 None"""
    settings = mock_ctx.config.settings
    settings.mode = "local"
    webhook_url = settings.integration.webhook.url

    # 注册 responses 拦截
    responses.add(responses.POST, webhook_url, status=200)

    send_webhook(mock_ctx)

    assert len(responses.calls) == 1
    request = responses.calls[0].request

    # 验证 JSON 载荷
    payload = json.loads(request.body)  # type: ignore
    assert payload["swanlab"]["mode"] == "local"
    assert payload["swanlab"]["exp_url"] is None
