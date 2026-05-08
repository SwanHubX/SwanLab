"""
Unit tests for LarkCallback and NotificationCallback base class.
Mocks HTTP layer — no real network calls.
"""

import base64
import hashlib
import hmac
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from swanlab.plugin.notification import NotificationCallback
from swanlab.plugin.notification.base import build_run_url
from swanlab.plugin.notification.lark import LarkCallback, _LarkBot

# ---------------------------------------------------------------------------
# Helpers & fixtures
# ---------------------------------------------------------------------------


def _mock_settings(**overrides):
    """Create a real Settings instance suitable for isinstance checks."""
    from swanlab.sdk.internal.settings import Settings

    s = Settings()
    overrides.setdefault("mode", "online")
    for k, v in overrides.items():
        if "." in k:
            parts = k.split(".", 1)
            sub = getattr(s, parts[0])
            object.__setattr__(sub, parts[1], v)
        else:
            object.__setattr__(s, k, v)
    return s


def _mock_response(json_body=None, raise_on_status=False):
    r = MagicMock()
    r.json.return_value = json_body or {"code": 0}
    r.raise_for_status.side_effect = requests.HTTPError("bad") if raise_on_status else None
    return r


@pytest.fixture
def lark_cb():
    cb = LarkCallback(webhook_url="https://lark.example.com/webhook")
    cb._settings = _mock_settings()
    cb._path = "/test-user/test-project/abc123"
    return cb


# ---------------------------------------------------------------------------
# _LarkBot
# ---------------------------------------------------------------------------


class TestLarkBot:
    def test_gen_sign(self):
        bot = _LarkBot("https://example.com/webhook", secret="my-secret")
        ts = 1700000000
        sign = bot.gen_sign(ts)
        expected = base64.b64encode(hmac.new(f"{ts}\nmy-secret".encode(), digestmod=hashlib.sha256).digest()).decode()
        assert sign == expected

    def test_gen_sign_without_secret_raises(self):
        bot = _LarkBot("https://example.com/webhook")
        with pytest.raises(ValueError, match="secret is required"):
            bot.gen_sign(1700000000)

    def test_attributes(self):
        bot = _LarkBot("https://example.com/webhook", secret="s")
        assert bot.webhook_url == "https://example.com/webhook"
        assert bot.secret == "s"


# ---------------------------------------------------------------------------
# build_run_url
# ---------------------------------------------------------------------------


class TestBuildRunUrl:
    def test_online_with_path(self):
        s = _mock_settings(mode="online", web_host="https://swanlab.cn")
        url = build_run_url(s, "/test-user/test-project/abc123")
        assert url is not None
        assert url.startswith("https://swanlab.cn/@")

    def test_offline_with_path(self):
        s = _mock_settings(mode="offline", web_host="https://swanlab.cn")
        url = build_run_url(s, "/test-user/test-project/abc123")
        assert url is not None

    def test_local_returns_none(self):
        s = _mock_settings(mode="local")
        assert build_run_url(s, "/test-user/test-project/abc123") is None

    def test_no_path_returns_none(self):
        s = _mock_settings(mode="online")
        assert build_run_url(s, None) is None


# ---------------------------------------------------------------------------
# NotificationCallback base
# ---------------------------------------------------------------------------


class TestNotificationCallbackBase:
    def test_name_returns_class_name(self):
        class MyCb(NotificationCallback):
            def _send_notification(self, state, error): ...

        assert MyCb().name == "MyCb"

    def test_on_run_initialized_captures_settings_and_path(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        settings = _mock_settings()
        cb.on_run_initialized(Path("/tmp"), "/test-user/test-project/abc123", settings=settings)
        assert cb._settings is settings
        assert cb._path == "/test-user/test-project/abc123"

    def test_on_run_initialized_ignores_non_settings(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb.on_run_initialized(Path("/tmp"), "/path", settings="not-a-Settings")
        assert cb._settings is None
        assert cb._path is None

    def test_on_run_initialized_no_kwarg(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb.on_run_initialized(Path("/tmp"), "/path")
        assert cb._settings is None
        assert cb._path is None

    def test_on_run_finished_noop_without_settings(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        with patch.object(cb._executor, "run") as mock_run:
            cb.on_run_finished(state="finished")
            mock_run.assert_not_called()

    def test_on_run_finished_submits_to_executor(self, lark_cb):
        with patch.object(lark_cb._executor, "run") as mock_run:
            lark_cb.on_run_finished(state="finished")
            mock_run.assert_called_once()
            assert mock_run.call_args[0][0] == lark_cb._send_notification
            assert mock_run.call_args[0][1] == "finished"

    def test_on_run_finished_passes_error(self, lark_cb):
        with patch.object(lark_cb._executor, "run") as mock_run:
            lark_cb.on_run_finished(state="crashed", error="OOM")
            assert mock_run.call_args[0][2] == "OOM"

    @patch("swanlab.plugin.notification.base.build_run_url", return_value="https://swanlab.cn/@u/p/runs/1")
    def test_get_url_with_path(self, _mock_build):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb._settings = _mock_settings()
        cb._path = "/u/p/1"
        assert cb._get_url() == "https://swanlab.cn/@u/p/runs/1"

    def test_get_url_without_settings(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        assert cb._get_url() is None


# ---------------------------------------------------------------------------
# LarkCallback._send_notification (mock HTTP)
# ---------------------------------------------------------------------------


class TestLarkCallbackSend:
    @patch("swanlab.plugin.notification.lark.requests.post")
    @patch("swanlab.plugin.notification.lark.console")
    def test_success_without_secret(self, mock_console, mock_post, lark_cb):
        mock_post.return_value = _mock_response()
        lark_cb._send_notification("finished", None)

        payload = mock_post.call_args[1]["json"]
        assert payload["msg_type"] == "text"
        assert "sign" not in payload
        assert mock_post.call_args[0][0] == "https://lark.example.com/webhook"
        mock_console.info.assert_called_once()
        mock_console.warning.assert_not_called()

    @patch("swanlab.plugin.notification.lark.requests.post")
    def test_success_with_secret_includes_sign(self, mock_post):
        cb = LarkCallback(webhook_url="https://lark.example.com/webhook", secret="my-secret")
        cb._settings = _mock_settings()
        cb._path = "/test-user/test-project/abc123"
        mock_post.return_value = _mock_response()
        cb._send_notification("finished", None)

        payload = mock_post.call_args[1]["json"]
        assert "sign" in payload
        assert payload["timestamp"] > 0

    @patch("swanlab.plugin.notification.lark.requests.post")
    @patch("swanlab.plugin.notification.lark.console")
    def test_business_error_warns(self, mock_console, mock_post, lark_cb):
        mock_post.return_value = _mock_response(json_body={"code": 19001, "msg": "invalid sign"})
        lark_cb._send_notification("finished", None)

        mock_console.warning.assert_called_once()
        assert "LarkBot sending failed" in mock_console.warning.call_args[0][0]
        assert "invalid sign" in mock_console.warning.call_args[0][0]

    @patch("swanlab.plugin.notification.lark.requests.post")
    def test_payload_contains_error_content(self, mock_post):
        cb = LarkCallback(webhook_url="https://lark.example.com/webhook", language="en")
        cb._settings = _mock_settings()
        cb._path = "/test-user/test-project/abc123"
        mock_post.return_value = _mock_response()
        cb._send_notification("crashed", "OOM")

        text = mock_post.call_args[1]["json"]["content"]["text"]
        assert "OOM" in text
        assert "test-project" in text
