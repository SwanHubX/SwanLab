"""
@author: cunyue
@file: test_login.py
@time: 2026/3/8 17:30
@description: 测试 SwanLab 登录流程 (E2E)
"""

from unittest.mock import patch

import pytest
import responses

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd.login import login
from swanlab.sdk.internal.context import RunConfig, RunContext, set_context
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.settings import settings


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

    @responses.activate
    def test_login_success_with_explicit_key(self):
        """测试：显式传入 API Key 并成功完成网络请求登录"""
        api_key = "test-explicit-key"

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

        # 这里的 patch 路径改为了具体的模块引用路径
        with patch("swanlab.sdk.cmd.login.apikey.prompt", return_value=prompt_key) as mock_prompt:
            result = login(api_key=None, save=True)

            # 验证全局单例被正确更新
            assert settings.web_host == "https://swanlab.cn"
            assert settings.api_host == "https://api.swanlab.cn"
            assert settings.api_key == prompt_key

            assert result is True
            mock_prompt.assert_called_once()

            # 验证凭证保存到了物理文件
            nrc_path = tmp_path / ".netrc"
            assert nrc_path.exists()
            content = nrc_path.read_text()
            assert prompt_key in content

    @responses.activate
    def test_login_custom_host_and_relogin(self):
        """测试：传入自定义 host，并且强制重新登录"""
        custom_host = "http://private-swanlab.local"
        custom_key = "private-key"

        responses.add(
            responses.POST,
            f"{custom_host}/api/login/api_key",
            json={"sid": "mock_private_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

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
        # 模拟 401 鉴权失败
        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"code": 401, "message": "Invalid API Key"},
            status=401,
        )

        # 修复：捕获你在 login.py 中抛出的自定义 AuthenticationError
        with pytest.raises(AuthenticationError, match="Failed to login"):
            login(api_key="wrong-key")
