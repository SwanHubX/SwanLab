"""
@author: cunyue
@file: test_login_e2e.py
@time: 2026/3/14
@description: 测试 SwanLab 登录流程 (E2E)
"""

import pytest
import responses

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd.login import login
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.settings import settings


class TestLoginE2E:
    @responses.activate
    def test_login_host_cleaning(self):
        """测试：传入不规范的 host 时，SDK 应自动清洗（补全协议、去除路径和 query）"""
        messy_host = "  10.0.0.1:8080/api/v1/?token=123  "
        expected_clean_host = "https://10.0.0.1:8080"
        api_key = "test-key"

        responses.add(
            responses.POST,
            f"{expected_clean_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        result = login(api_key=api_key, host=messy_host)

        assert result is True
        assert settings.api_host == expected_clean_host
        assert settings.web_host == expected_clean_host

    @responses.activate
    def test_login_official_host_protection(self):
        """测试：传入官方 api_host 时，不应错误覆盖 web_host"""
        official_host = "api.swanlab.cn"
        expected_api_host = "https://api.swanlab.cn"
        api_key = "test-key"

        responses.add(
            responses.POST,
            f"{expected_api_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        assert settings.web_host == "https://swanlab.cn"

        result = login(api_key=api_key, host=official_host, relogin=True)

        assert result is True
        assert settings.api_host == expected_api_host
        # 官方 web_host 不得被 api 子域名覆盖
        assert settings.web_host == "https://swanlab.cn"

    def test_login_skip_if_already_logged_in(self, monkeypatch):
        """测试：如果已经登录且没有强制 relogin，应直接跳过并返回 True"""

        monkeypatch.setattr("swanlab.sdk.cmd.login.client.exists", lambda: True)
        result = login(api_key="some_key", relogin=False)
        assert result is True

    def test_login_block_if_run_active(self):
        """测试：如果 SwanLab Run 正在运行，不允许登录"""
        from swanlab.sdk.cmd.init import init
        from swanlab.sdk.internal.run import has_run

        # 先创建一个真实的 run
        init(mode="disabled")
        assert has_run()

        # 尝试登录应该抛出异常
        with pytest.raises(RuntimeError, match="`swanlab.login` requires no active SwanLabRun"):
            login(api_key="some_key")

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
    def test_login_with_prompt_and_save(self, monkeypatch, tmp_path):
        """测试：没有 API Key 时触发交互式输入，并测试凭证保存 (save=True)"""
        prompt_key = "test-prompt-key"
        prompt_called = []

        responses.add(
            responses.POST,
            "https://fake.swanlab.cn/api/login/api_key",
            json={"sid": "mock_sid", "user": "mock_user", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        monkeypatch.setattr("swanlab.sdk.cmd.login.apikey.prompt", lambda **_: (prompt_called.append(1), prompt_key)[1])

        result = login(api_key=None, host="fake.swanlab.cn", save=True)

        assert settings.web_host == "https://fake.swanlab.cn"
        assert settings.api_host == "https://fake.swanlab.cn"
        assert settings.api_key == prompt_key
        assert result is True
        assert len(prompt_called) == 1

        nrc_path = settings.root / ".netrc"
        assert nrc_path.exists()
        content = nrc_path.read_text()
        assert prompt_key in content

    @responses.activate
    def test_login_custom_host_and_relogin(self, monkeypatch):
        """测试：传入自定义 host，并且强制重新登录"""
        custom_host = "http://private-swanlab.local"
        custom_key = "private-key"

        responses.add(
            responses.POST,
            f"{custom_host}/api/login/api_key",
            json={"sid": "mock_private_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        monkeypatch.setattr("swanlab.sdk.cmd.login.client.exists", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.login.client.reset", lambda: None)

        result = login(api_key=custom_key, host=custom_host, relogin=True)

        assert result is True
        assert len(responses.calls) == 1
        assert responses.calls[0].request.url == f"{custom_host}/api/login/api_key"
        assert settings.api_host == custom_host
        assert settings.api_key == custom_key

    @responses.activate
    def test_login_network_failure(self):
        """测试：网络请求失败或者 API Key 错误时抛出 AuthenticationError"""
        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"code": 401, "message": "Invalid API Key"},
            status=401,
        )

        with pytest.raises(AuthenticationError, match="Failed to login"):
            login(api_key="wrong-key", relogin=True)

    @responses.activate
    def test_login_host_changed_triggers_prompt(self, monkeypatch):
        """测试：上次保存了旧 host 的 API Key，本次传入新 host 时不复用旧 key，应触发重新输入"""
        old_host = "https://old.swanlab.cn"
        old_prompt_key = "old-private-key"
        new_host = "https://private.swanlab.com"
        new_prompt_key = "new-private-key"

        responses.add(
            responses.POST,
            f"{old_host}/api/login/api_key",
            json={"sid": "mock_old_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )
        responses.add(
            responses.POST,
            f"{new_host}/api/login/api_key",
            json={"sid": "mock_new_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        prompt_calls = []

        def mock_prompt(**_):
            if len(prompt_calls) == 0:
                prompt_calls.append(old_prompt_key)
                return old_prompt_key
            else:
                prompt_calls.append(new_prompt_key)
                return new_prompt_key

        monkeypatch.setattr("swanlab.sdk.cmd.login.apikey.prompt", mock_prompt)

        login(api_key=None, host=old_host, save=True)
        assert len(responses.calls) == 1

        result = login(api_key=None, host=new_host, relogin=True)

        assert result is True
        assert len(prompt_calls) == 2
        assert len(responses.calls) == 2
        assert responses.calls[1].request.headers["authorization"] == new_prompt_key
