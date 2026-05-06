"""
Unit tests for LarkCallback and _NotificationCallback base class.
Mocks HTTP layer — no real network calls.
"""

import base64
import hashlib
import hmac
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from swanlab.plugin.notification import _NotificationCallback
from swanlab.plugin.notification.lark import LarkCallback, _LarkBot
from swanlab.sdk.typings.run.callback import RunInfo

# ---------------------------------------------------------------------------
# Helpers & fixtures
# ---------------------------------------------------------------------------


def _make_run_info(**overrides) -> RunInfo:
    defaults = dict(
        project="test-project",
        workspace="test-user",
        experiment_name="exp-001",
        description="a test experiment",
        run_id="abc123",
        run_dir=Path("/tmp/swanlab"),
        path="/test-user/test-project/abc123",
        mode="online",
        url="https://swanlab.cn/@test-user/test-project/runs/abc123",
    )
    defaults.update(overrides)
    return RunInfo(**defaults)  # type: ignore


def _mock_response(json_body=None, raise_on_status=False):
    r = MagicMock()
    r.json.return_value = json_body or {"code": 0}
    r.raise_for_status.side_effect = requests.HTTPError("bad") if raise_on_status else None
    return r


@pytest.fixture
def lark_cb():
    cb = LarkCallback(webhook_url="https://lark.example.com/webhook")
    cb._run_info = _make_run_info()
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
# _NotificationCallback base
# ---------------------------------------------------------------------------


class TestNotificationCallbackBase:
    def test_name_returns_class_name(self):
        class MyCb(_NotificationCallback):
            def _send_notification(self, state, error): ...

        assert MyCb().name == "MyCb"

    def test_on_run_initialized_extracts_run_info(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        info = _make_run_info()
        cb.on_run_initialized(Path("/tmp"), "/test-user/test-project/abc123", run_info=info)
        assert cb._run_info is info

    def test_on_run_initialized_ignores_non_run_info(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb.on_run_initialized(Path("/tmp"), "/path", run_info="not-a-RunInfo")
        assert cb._run_info is None

    def test_on_run_initialized_no_kwarg(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb.on_run_initialized(Path("/tmp"), "/path")
        assert cb._run_info is None

    def test_on_run_finished_noop_without_run_info(self):
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

    def test_get_url_with_run_info(self):
        cb = LarkCallback(webhook_url="https://example.com/webhook")
        cb._run_info = _make_run_info(url="https://swanlab.cn/@u/p/runs/1")
        assert cb._get_url() == "https://swanlab.cn/@u/p/runs/1"

    def test_get_url_without_run_info(self):
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
        # posts to correct webhook url
        assert mock_post.call_args[0][0] == "https://lark.example.com/webhook"
        mock_console.info.assert_called_once()
        mock_console.warning.assert_not_called()

    @patch("swanlab.plugin.notification.lark.requests.post")
    def test_success_with_secret_includes_sign(self, mock_post):
        cb = LarkCallback(webhook_url="https://lark.example.com/webhook", secret="my-secret")
        cb._run_info = _make_run_info()
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
    def test_http_error_raises(self, mock_post, lark_cb):
        mock_post.side_effect = requests.ConnectionError("timeout")

        with pytest.raises(requests.ConnectionError):
            lark_cb._send_notification("finished", None)

    @patch("swanlab.plugin.notification.lark.requests.post")
    def test_payload_contains_error_content(self, mock_post):
        cb = LarkCallback(webhook_url="https://lark.example.com/webhook", language="en")
        cb._run_info = _make_run_info()
        mock_post.return_value = _mock_response()

        cb._send_notification("crashed", "OOM")

        text = mock_post.call_args[1]["json"]["content"]["text"]
        assert "OOM" in text
        assert "test-project" in text
