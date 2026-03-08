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
    @responses.activate
    def test_login_host_cleaning(self):
        """测试：传入不规范的 host 时，SDK 应自动清洗（补全协议、去除路径和 query）"""
        # 带有前后空格、没有 http 头、带有 path 和斜杠的“脏数据”
        messy_host = "  10.0.0.1:8080/api/v1/?token=123  "
        expected_clean_host = "https://10.0.0.1:8080"
        api_key = "test-key"

        # 验证清洗后的 URL 是否被正确用于请求
        responses.add(
            responses.POST,
            f"{expected_clean_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 这里用 relogin=True 绕过已登录状态检查
        result = login(api_key=api_key, host=messy_host)

        assert result is True
        # 验证 api_host 被清洗干净
        assert settings.api_host == expected_clean_host
        # 验证私有化部署下，web_host 跟随 api_host 更新
        assert settings.web_host == expected_clean_host

    @responses.activate
    def test_login_official_host_protection(self):
        """测试：传入官方 api_host 时，不应错误覆盖 web_host"""
        # 用户可能只传入了 api 子域名
        official_host = "api.swanlab.cn"
        expected_api_host = "https://api.swanlab.cn"
        api_key = "test-key"

        responses.add(
            responses.POST,
            f"{expected_api_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 假装当前环境未被污染，默认 web_host 是正确的
        assert settings.web_host == "https://swanlab.cn"

        result = login(api_key=api_key, host=official_host, relogin=True)

        assert result is True
        # 验证 api_host 正常更新了协议头
        assert settings.api_host == expected_api_host
        # 【核心断言】：验证官方 web_host 没有被改成 api 域名，防线生效
        assert settings.web_host == "https://swanlab.cn"

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
            with patch("swanlab.sdk.cmd.login.client.reset"):
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
