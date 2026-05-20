"""
@author: cunyue
@file: test_sync_helper.py.py
@time: 2026/5/17 19:30
@description: 测试sync的工具函数
"""

import pytest
from pydantic import ValidationError

from swanlab.exceptions import AuthenticationError
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncFlushResponse, DeliverSyncStartResponse
from swanlab.sdk.cmd.sync import ensure_run_dir, sync
from swanlab.sdk.internal.settings import Settings


def test_ensure_run_dir_ok(tmp_path):
    assert ensure_run_dir(tmp_path) == tmp_path.resolve()


def test_ensure_run_dir_requires_directory(tmp_path):
    file_path = tmp_path / "run.txt"
    file_path.write_text("", encoding="utf-8")

    with pytest.raises(ValidationError):
        ensure_run_dir(file_path)


def test_ensure_run_dir_requires_readable(tmp_path, monkeypatch):
    monkeypatch.setattr("swanlab.sdk.cmd.sync.os.access", lambda *_: False)

    with pytest.raises(PermissionError):
        ensure_run_dir(tmp_path)


def test_sync_requires_api_key_when_client_missing(tmp_path, monkeypatch):
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: False)

    with pytest.raises(AuthenticationError, match="Please login first"):
        sync(tmp_path, settings=Settings())


def test_sync_logs_in_when_api_key_provided(tmp_path, monkeypatch):
    logged_in = False
    login_calls = []

    def client_exists():
        return logged_in

    def login_raw(api_key, *, host, print_welcome):
        nonlocal logged_in
        login_calls.append((api_key, host, print_welcome))
        logged_in = True

    class FakeCoreSync:
        def confirm_sync_finish(self):
            return type("Response", (), {"success": True, "message": "OK"})()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", client_exists)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.login_raw", login_raw)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: FakeCoreSync())
    monkeypatch.setattr("swanlab.sdk.cmd.sync._deliver_sync_start", lambda *args, **kwargs: None)

    sync(tmp_path, settings=Settings(api_key="test-api-key", api_host="https://api.example.com"))

    assert login_calls == [("test-api-key", "https://api.example.com", False)]


def test_sync_raises_when_deliver_sync_start_fails(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=False, message="start failed")

        def deliver_sync_flush(self):
            raise AssertionError("flush should not be called")

        def confirm_sync_finish(self):
            self.confirmed = True
            return type("Response", (), {"success": True, "message": "OK"})()

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="start failed"):
        sync(tmp_path, settings=Settings())

    assert core.confirmed is False


def test_sync_raises_when_deliver_sync_flush_fails(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=True, message="success")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=False, message="flush failed")

        def confirm_sync_finish(self):
            self.confirmed = True
            return type("Response", (), {"success": True, "message": "OK"})()

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="flush failed"):
        sync(tmp_path, settings=Settings())

    assert core.confirmed is False
