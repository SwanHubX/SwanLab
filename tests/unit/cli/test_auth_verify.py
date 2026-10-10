import pytest
import responses
from click.testing import CliRunner

from swanlab.cli.auth.verify import verify
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console, nrc
from swanlab.sdk.internal.settings import Settings


@pytest.fixture()
def saved_credentials():
    path = Settings.get_user_config_dir() / ".netrc"
    path.parent.mkdir(parents=True, exist_ok=True)
    nrc.write(path, api_host="https://server.example", web_host="https://server.example", api_key="test-key")


@pytest.mark.parametrize(
    ("status", "message"),
    [
        (403, "user is not verified"),
        (404, "Please upgrade"),
        (429, "rate limit exceeded"),
        (503, "service unavailable"),
    ],
)
@responses.activate
def test_verify_reports_error_category(status, message, monkeypatch, saved_credentials):
    responses.add(responses.GET, "https://server.example/api/auth/verify", status=status, json={"message": message})
    errors = []
    monkeypatch.setattr(console, "error", errors.append)
    result = CliRunner().invoke(verify)
    assert result.exit_code == 1
    assert message in errors[0]
    assert not client.exists()
    assert len(responses.calls) == 1


@responses.activate
def test_verify_success_exits_zero(monkeypatch, saved_credentials):
    responses.add(
        responses.GET,
        "https://server.example/api/auth/verify",
        json={"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"},
    )
    errors = []
    messages = []
    monkeypatch.setattr(console, "error", errors.append)
    monkeypatch.setattr(console, "info", lambda *args, **kwargs: messages.append(args))
    result = CliRunner().invoke(verify)
    assert result.exit_code == 0
    assert not errors
    assert str(messages[0][-1]) == "alice"
    assert not client.exists()
