from threading import Event, Thread
from unittest.mock import Mock

import pytest

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, RunState
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.core import CorePython
from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.core_python.transport.thread import Transport
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


@pytest.fixture(params=["run", "sync"])
def owned_core(request, core_auth, monkeypatch, tmp_path):
    core = CorePython("online") if request.param == "run" else CoreSyncPython()
    core._ctx = CoreContext(
        config=CoreConfig(
            run_id="lifecycle",
            run_dir=tmp_path,
            section_rule=0,
            record_batch=10,
            record_interval=0.01,
            save_split=1024,
            save_size=1024,
            save_part=1024,
            save_batch=10,
        )
    )
    core._ctx.set_online_params("test-user", "demo", "project", 1, "experiment")
    core._client_instance = client.new("core-test-key", core_auth)
    core._tracker = UploadTracker()
    core._tracker.set_state(CoreState.CORE_STATE_RUNNING)
    if isinstance(core, CoreSyncPython):
        monkeypatch.setattr(core, "_read_executor", Mock())
        monkeypatch.setattr(core, "_reader", Mock())
    yield core
    if isinstance(core, CorePython):
        core._rollback_run_start()
    else:
        core._rollback_sync_start()


@pytest.mark.parametrize(
    ("interrupted", "stop_error"),
    [
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    ],
    ids=["timeout", "stop_error", "interrupted", "interrupted_stop_error"],
)
def test_active_worker_retains_client_until_cleanup_succeeds(
    owned_core, core_auth, monkeypatch, interrupted, stop_error
):
    core = owned_core
    ctx = core._ctx
    owned = core._client_instance
    entered, release = Event(), Event()
    observed_clients = []

    def dispatch(_records):
        entered.set()
        release.wait()
        observed_clients.append(client._get_client())
        return True, []

    transport = Transport(ctx, auto_start=False)
    monkeypatch.setattr(transport, "_dispatcher", dispatch)
    transport._thread = Thread(target=transport._loop, daemon=True)
    transport._thread.start()
    core._transport = transport
    monkeypatch.setattr(Transport, "FINISH_JOIN_TIMEOUT", 0.01)
    finish = transport.finish

    def fail_stop(timeout=None):
        if interrupted and timeout is None:
            raise KeyboardInterrupt()
        if stop_error:
            raise RuntimeError("join failed")
        return finish(timeout)

    transport.put([Record()])
    assert entered.wait(2)
    try:
        monkeypatch.setattr(transport, "finish", fail_stop)
        if interrupted:
            if isinstance(core, CorePython):
                action = core._confirm_finish_when_enabled
            else:
                core._read_executor.wait.side_effect = KeyboardInterrupt()
                action = core.confirm_sync_finish
            with pytest.raises(KeyboardInterrupt):
                action()
        elif isinstance(core, CorePython):
            core._rollback_run_start()
        else:
            core._rollback_sync_start()
        assert transport.is_alive()
        assert core._transport is transport
        assert core._tracker is not None
        assert core._tracker.snapshot().state == CoreState.CORE_STATE_RUNNING
        assert core._client_instance is not None
        assert client._get_client() is owned
        with pytest.raises(RuntimeError, match="already exists"):
            client.new("replacement", core_auth)
    finally:
        monkeypatch.setattr(transport, "finish", finish)
        release.set()
        assert transport.finish(timeout=2)
        if isinstance(core, CorePython):
            core._rollback_run_start()
        else:
            core._rollback_sync_start()
    assert observed_clients == [owned]
    assert not client.exists()
    assert core._client_instance is None


def test_other_consumer_cleanup_error_retains_ownership(owned_core, monkeypatch):
    core = owned_core
    owned = core._client_instance
    consumer = Mock()
    if isinstance(core, CorePython):
        monkeypatch.setattr(core, "_heartbeat", consumer)
        consumer.stop.side_effect = RuntimeError("heartbeat still active")
        rollback = core._rollback_run_start
    else:
        monkeypatch.setattr(core, "_read_executor", consumer)
        monkeypatch.setattr(core, "_reader", Mock())
        consumer.close.side_effect = RuntimeError("reader still active")
        rollback = core._rollback_sync_start
    try:
        rollback()
        assert client._get_client() is owned
        assert core._client_instance is not None
        if isinstance(core, CorePython):
            assert core._heartbeat is consumer
    finally:
        consumer.stop.side_effect = None
        consumer.close.side_effect = None
        rollback()
    assert not client.exists()


def test_finish_reports_client_close_failure_and_can_retry(owned_core, monkeypatch):
    core = owned_core
    owned = core._client_instance
    transport = Mock()
    transport.finish.return_value = True
    core._transport = transport
    finish_record = FinishRecord(state=RunState.RUN_STATE_FINISHED)
    if isinstance(core, CorePython):
        core._pending_online_finish_record = finish_record
        confirm = core.confirm_run_finish
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", Mock())
    else:
        core._finish_record = finish_record
        confirm = core.confirm_sync_finish
        monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", Mock())
    with monkeypatch.context() as patcher:
        patcher.setattr(owned, "close", Mock(side_effect=RuntimeError("close failed")))
        response = confirm()
        assert not response.success
        assert "Failed to release" in response.message
        assert core._client_instance is owned
        assert client._get_client() is owned
    assert core._teardown_client()
    assert not client.exists()
