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
from swanlab.sdk.internal.pkg.netrc import write_netrc
from swanlab.sdk.internal.settings import Settings, settings


class TestLoginE2E:
    @staticmethod
    def _reinit_settings():
        """重新从环境变量、.netrc 等来源初始化全局 settings，模拟新会话启动。"""
        fresh = Settings()
        for field_name in Settings.model_fields.keys():
            object.__setattr__(settings, field_name, getattr(fresh, field_name))
        object.__setattr__(settings, "__pydantic_fields_set__", set())

    @staticmethod
    def _save_netrc(host: str, key: str):
        """将凭证写入 .netrc 并重新初始化 settings，模拟上次会话保存、新会话加载的完整流程。"""
        nrc_path = settings.root / ".netrc"
        nrc_path.parent.mkdir(parents=True, exist_ok=True)
        write_netrc(nrc_path, host=host, username=host, password=key)
        TestLoginE2E._reinit_settings()

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

    @responses.activate
    def test_login_skip_if_already_logged_in(self):
        """测试：如果已经登录且没有强制 relogin，无论是否传入入参，应直接跳过并返回 True"""
        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )
        login(api_key="first-key")
        assert client.exists()

        # 第二次调用不带 relogin=True，应直接跳过，不再发起网络请求
        result = login(relogin=False)
        assert result is True
        assert len(responses.calls) == 1

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
    def test_login_with_prompt_and_save(self, monkeypatch):
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
    def test_login_custom_host_and_relogin(self):
        """测试：已登录状态下，传入自定义 host 并强制重新登录"""
        initial_key = "initial-key"
        custom_host = "http://private-swanlab.local"
        custom_key = "private-key"

        responses.add(
            responses.POST,
            "https://api.swanlab.cn/api/login/api_key",
            json={"sid": "initial_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )
        responses.add(
            responses.POST,
            f"{custom_host}/api/login/api_key",
            json={"sid": "mock_private_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 先做一次真实登录，使 client 进入已登录状态
        login(api_key=initial_key)
        assert client.exists()

        result = login(api_key=custom_key, host=custom_host, relogin=True)

        assert result is True
        assert len(responses.calls) == 2
        assert responses.calls[1].request.url == f"{custom_host}/api/login/api_key"
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

    @responses.activate
    def test_repeated_login_with_different_credentials(self):
        """测试：同一会话中先后调用 login，每次传入不同的 api_key 和 host，应以最新入参为准重新登录"""
        host_a = "https://server-a.swanlab.com"
        host_b = "https://server-b.swanlab.com"
        key_a = "key-for-server-a"
        key_b = "key-for-server-b"

        responses.add(
            responses.POST,
            f"{host_a}/api/login/api_key",
            json={"sid": "sid_a", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )
        responses.add(
            responses.POST,
            f"{host_b}/api/login/api_key",
            json={"sid": "sid_b", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 第一次登录
        result = login(api_key=key_a, host=host_a)
        assert result is True
        assert settings.api_host == host_a
        assert settings.api_key == key_a

        # 第二次登录，传入不同凭证
        result = login(api_key=key_b, host=host_b, relogin=True)
        assert result is True
        assert len(responses.calls) == 2
        assert responses.calls[1].request.headers["authorization"] == key_b
        assert settings.api_host == host_b
        assert settings.api_key == key_b

    @responses.activate
    def test_login_explicit_params_override_netrc(self):
        """测试：本地 .netrc 中已保存旧凭证（模拟上次会话），新会话首次登录时显式传入不同的
        api_key 和 host，应优先使用入参，而非 .netrc 中的旧凭证"""
        old_host = "https://old.swanlab.com"
        old_key = "old-key-from-netrc"
        new_host = "https://new.swanlab.com"
        new_key = "new-key"

        # 写入旧凭证到 .netrc（settings.root 由 conftest 指向 tmp_path/.swanlab）
        self._save_netrc(old_host, old_key)

        # 前提：settings 已从 .netrc 加载旧凭证
        assert settings.api_key == old_key
        assert settings.api_host == old_host

        responses.add(
            responses.POST,
            f"{new_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 传入新凭证登录 —— 入参应覆盖 .netrc 中加载的旧凭证
        result = login(api_key=new_key, host=new_host)

        assert result is True
        assert len(responses.calls) == 1
        assert responses.calls[0].request.headers["authorization"] == new_key
        assert settings.api_host == new_host
        assert settings.api_key == new_key

    @responses.activate
    def test_login_autouse_netrc_credentials(self):
        """测试：本地 .netrc 中已保存凭证，login() 不传任何入参时，应自动使用存储的 key 完成登录，
        不弹 prompt"""
        stored_host = "https://stored.swanlab.com"
        stored_key = "stored-key-from-netrc"

        # 写入凭证到 .netrc 并重新初始化 settings，模拟新会话从 .netrc 自动加载
        self._save_netrc(stored_host, stored_key)

        assert settings.api_key == stored_key
        assert settings.api_host == stored_host

        responses.add(
            responses.POST,
            f"{stored_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 不传任何入参，应自动使用 .netrc 中的凭证
        result = login()

        assert result is True
        assert len(responses.calls) == 1
        assert responses.calls[0].request.headers["authorization"] == stored_key
        assert settings.api_host == stored_host
        assert settings.api_key == stored_key

    @responses.activate
    def test_login_env_key_reused_when_host_changes_without_netrc(self, monkeypatch):
        """测试：环境变量中有 api_key，login 只传新 host 不传 api_key，且本地无 .netrc 时，
        应静默复用环境变量中的 key（不弹 prompt，不打 warning）"""
        env_key = "key-from-env"
        env_host = "https://env.swanlab.com"
        new_host = "https://new.swanlab.com"

        # 设置环境变量，再重建 Settings 实例模拟程序启动时的自动解析
        monkeypatch.setenv("SWANLAB_API_KEY", env_key)
        monkeypatch.setenv("SWANLAB_API_HOST", env_host)
        self._reinit_settings()

        # 前提：env var 已加载，本地无 .netrc 文件
        assert settings.api_key == env_key
        assert settings.api_host == env_host
        assert not (settings.root / ".netrc").exists()

        responses.add(
            responses.POST,
            f"{new_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        # 只传新 host，不传 api_key —— 无 .netrc，env key 应被静默复用
        result = login(host=new_host)

        assert result is True
        assert len(responses.calls) == 1
        assert responses.calls[0].request.headers["authorization"] == env_key
        assert settings.api_host == new_host

    @responses.activate
    def test_login_explicit_params_override_env_vars(self, monkeypatch):
        """测试：环境变量中已配置 SWANLAB_API_KEY / SWANLAB_API_HOST（settings 自动解析），
        调用 login 时传入不同的 api_key 和 host，应优先使用入参，而非环境变量中的值"""
        env_key = "key-from-env"
        env_host = "https://env.swanlab.com"
        new_host = "https://custom.swanlab.com"
        new_key = "custom-key"

        monkeypatch.setenv("SWANLAB_API_KEY", env_key)
        monkeypatch.setenv("SWANLAB_API_HOST", env_host)
        self._reinit_settings()

        assert settings.api_key == env_key
        assert settings.api_host == env_host

        responses.add(
            responses.POST,
            f"{new_host}/api/login/api_key",
            json={"sid": "mock_sid", "expiredAt": "2099-12-31T23:59:59.000Z"},
            status=200,
        )

        result = login(api_key=new_key, host=new_host)

        assert result is True
        assert len(responses.calls) == 1
        assert responses.calls[0].request.headers["authorization"] == new_key
        assert settings.api_host == new_host
        assert settings.api_key == new_key
