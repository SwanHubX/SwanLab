"""
@author: cunyue
@file: test_login.py
@time: 2026/3/8 17:30
@description: 测试 SwanLab 登录流程 (E2E)
"""

from unittest.mock import patch

import pytest
import responses

from swanlab.sdk.cmd.login import login
from swanlab.sdk.internal.context import RunConfig, RunContext, clear_context, set_context
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.settings import settings


@pytest.fixture(autouse=True)
def cleanup_swanlab(monkeypatch, tmp_path):
    # 1. 路径隔离（最重要）：利用 pytest 自带的 tmp_path，配合环境变量隔离
    # monkeypatch 的好处是测试一结束它会自动恢复原样，不需要手动 teardown
    monkeypatch.setenv("SWANLAB_SAVE_DIR", str(tmp_path))
    monkeypatch.delenv("SWANLAB_API_KEY", raising=False)
    monkeypatch.delenv("SWANLAB_API_HOST", raising=False)
    monkeypatch.delenv("SWANLAB_WEB_HOST", raising=False)

    yield  # 这里是执行测试用例的地方

    # --- 下面是手动 Teardown 阶段 ---

    # 2. 释放全局 Client（防止下一个测试报 "Client already exists"）
    from swanlab.sdk.internal.core_python import client

    if client.exists():
        client.reset()

    # 3. 清理运行上下文（防止 use_temp_context 报错或残留）
    from swanlab.sdk.internal.context import clear_context

    clear_context()


@pytest.fixture(autouse=True)
def isolate_env(tmp_path, monkeypatch):
    """
    全局环境隔离固件：
    1. 将 SWANLAB_SAVE_DIR 指向临时目录，隔离真实的 .netrc 和本地配置。
    2. 确保每次测试前，全局的 client 状态被彻底清空。
    """
    monkeypatch.setenv("SWANLAB_SAVE_DIR", str(tmp_path))
    monkeypatch.delenv("SWANLAB_API_KEY", raising=False)
    monkeypatch.delenv("SWANLAB_API_HOST", raising=False)

    # 测试前如果存在残留，清理掉
    if client.exists():
        client.reset()

    yield

    # 测试结束后的清理工作
    if client.exists():
        client.reset()
    clear_context()


class TestLoginE2E:
    def test_login_skip_if_already_logged_in(self):
        """测试：如果已经登录且没有强制 relogin，应直接跳过并返回 True"""
        with patch("swanlab.sdk.cmd.login.client.exists", return_value=True):
            result = login(api_key="some_key", relogin=False)
            assert result is True

    def test_login_block_if_context_active(self):
        """测试：如果运行上下文已存在（比如实验进行中），不允许登录"""
        mock_config = RunConfig(settings=settings)
        set_context(RunContext(config=mock_config))

        result = login(api_key="some_key")
        assert result is False

        clear_context()

    @responses.activate
    def test_login_success_with_explicit_key(self):
        """测试：显式传入 API Key 并成功完成网络请求登录"""
        api_key = "test-explicit-key"

        # 注意：此处加上了 /api 前缀，并补充了合法的 expiredAt
        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"sid": "mock_sid", "user": "mock_user", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        result = login(api_key=api_key, save=False)

        assert result is True
        assert settings.api_key == api_key
        assert len(responses.calls) == 1
        assert responses.calls[0].request.headers["authorization"] == api_key
        assert client.exists()

    @responses.activate
    def test_login_with_prompt_and_save(self, tmp_path):
        """测试：没有 API Key 时触发交互式输入，并测试凭证保存 (save=True)"""
        prompt_key = "test-prompt-key"

        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"sid": "mock_sid", "user": "mock_user", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        with patch("swanlab.sdk.cmd.login.apikey.prompt", return_value=prompt_key) as mock_prompt:
            result = login(api_key=None, save=True)

            assert result is True
            mock_prompt.assert_called_once()

            # 验证凭证保存
            nrc_path = tmp_path / ".netrc"
            assert nrc_path.exists()
            content = nrc_path.read_text()
            assert prompt_key in content

    @responses.activate
    def test_login_custom_host_and_relogin(self):
        """测试：传入自定义 host，并且强制重新登录"""
        custom_host = "http://private-swanlab.local"
        custom_key = "private-key"

        # 拦截私有化节点的请求（加上 /api）
        responses.add(
            responses.POST,
            f"{custom_host}/api/login/api_key",
            json={"sid": "mock_private_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 为了测试 relogin=True，先造一个假象让系统以为已经登录了
        with patch("swanlab.sdk.cmd.login.client.exists", return_value=True):
            result = login(api_key=custom_key, host=custom_host, relogin=True)

            assert result is True
            assert len(responses.calls) == 1
            assert responses.calls[0].request.url == f"{custom_host}/api/login/api_key"

            # 验证配置的同步更新
            assert settings.api_host == custom_host
            assert settings.api_key == custom_key

    @responses.activate
    def test_login_network_failure(self):
        """测试：网络请求失败或者 API Key 错误的情况"""
        # 模拟 401
        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"code": 401, "message": "Invalid API Key"},
            status=401,
        )

        # 鉴权失败时 login_resp 为 None，login() 里面有一句 assert login_resp is not None
        # 所以这里应当捕获到 AssertionError
        with pytest.raises(AssertionError, match="Failed to get login response"):
            login(api_key="wrong-key")
