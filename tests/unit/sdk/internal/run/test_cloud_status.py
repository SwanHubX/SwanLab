"""Run 动态云端路径（path / url）及 finish 收尾同步测试。"""

import threading
from types import SimpleNamespace
from typing import Any, cast
from unittest.mock import MagicMock

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import (
    CLOUD_STATE_ATTACHED,
    CLOUD_STATE_FAILED,
    CLOUD_STATE_PENDING,
    CLOUD_STATE_UNSPECIFIED,
    CloudStatus,
    CoreState,
    OperationStats,
    UploadStatus,
)
from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import Run, greeting


def _make_online_run(core: MagicMock, *, path=None, mode: str = "online") -> Run:
    run = Run.__new__(Run)
    settings = SimpleNamespace(
        mode=mode,
        web_host="https://swanlab.cn",
        experiment=SimpleNamespace(name="test-exp"),
    )
    ctx = SimpleNamespace(
        config=SimpleNamespace(settings=settings, run_dir="/tmp/swanlab-test-run"),
        run_dir="/tmp/swanlab-test-run",
        run_file="/tmp/swanlab-test-run/run-test.swanlab",
    )
    run._ctx = cast(Any, ctx)
    run._core = core
    run._path = path
    run._state = "running"
    run._init_pid = fork.current_pid()
    run._components = MagicMock()
    run._probe = MagicMock()
    run._callbacker = MagicMock()
    run._fork_registered = False
    run._atexit_registered = False
    run._excepthook_installed = False
    run._sigint_installed = False
    run._api_lock = threading.RLock()
    return run


class TestDynamicCloudProperties:
    """Run.path / Run.url 动态轮询云端状态的行为。"""

    def test_pending_then_attached_syncs_and_caches(self):
        core = MagicMock()
        core.get_operation_stats.side_effect = [
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    cloud=CloudStatus(state=CLOUD_STATE_PENDING),
                ),
            ),
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    cloud=CloudStatus(state=CLOUD_STATE_ATTACHED, path="/alice/proj/run-1"),
                ),
            ),
        ]
        run = _make_online_run(core, path=None)

        with pytest.raises(ValueError, match="has not connected to SwanLab online service yet"):
            _ = run.path

        assert run.path == "/alice/proj/run-1"
        assert run.url == "https://swanlab.cn/@alice/proj/runs/run-1"
        # _path 已缓存后不再重复查询 Core
        assert run.path == "/alice/proj/run-1"
        assert core.get_operation_stats.call_count == 2

    def test_failed_cloud_state_raises_reason(self):
        core = MagicMock()
        core.get_operation_stats.return_value = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(
                state=CoreState.CORE_STATE_RUNNING,
                cloud=CloudStatus(state=CLOUD_STATE_FAILED, error="401 invalid token"),
            ),
        )
        run = _make_online_run(core, path=None)

        with pytest.raises(ValueError, match="401 invalid token"):
            _ = run.url

    def test_non_online_mode_raises_without_polling(self):
        core = MagicMock()
        run = _make_online_run(core, path=None, mode="offline")

        with pytest.raises(ValueError, match="not using SwanLab online mode"):
            _ = run.path
        core.get_operation_stats.assert_not_called()


class TestFinishCloudLinkSync:
    """Run.finish() 期间延迟拿到云端路径时在收尾后补打链接。"""

    def test_prints_cloud_link_when_attached_during_finish(self, monkeypatch, tmp_path):
        core = MagicMock()
        core.confirm_after_deliver_failure = True
        core.last_error = None
        core.deliver_run_finish.return_value = DeliverRunFinishResponse(success=True, message="OK")
        core.confirm_run_finish.return_value = ConfirmRunFinishResponse(success=True, message="OK")
        # 第 1 次：goodbye 前 sync；第 2 次：run_with_progress 初始检查（进入轮询）；第 3、4 次：轮询退出及收尾同步拿到 path
        core.get_operation_stats.side_effect = [
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    cloud=CloudStatus(state=CLOUD_STATE_PENDING),
                ),
            ),
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    cloud=CloudStatus(state=CLOUD_STATE_PENDING),
                ),
            ),
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_FINISHED,
                    cloud=CloudStatus(state=CLOUD_STATE_ATTACHED, path="/alice/proj/run-9"),
                ),
            ),
            GetOperationStatsResponse(
                success=True,
                stats=OperationStats(
                    state=CoreState.CORE_STATE_FINISHED,
                    cloud=CloudStatus(state=CLOUD_STATE_ATTACHED, path="/alice/proj/run-9"),
                ),
            ),
        ]

        run = _make_online_run(core, path=None)
        ctx = cast(Any, run._ctx)
        ctx.run_dir = tmp_path
        ctx.config.settings = SimpleNamespace(
            mode="online",
            web_host="https://swanlab.cn",
            core=SimpleNamespace(skip_store=False),
            run=SimpleNamespace(dir=tmp_path.name),
            experiment=SimpleNamespace(name="test-exp"),
        )
        ctx.config.run_dir = tmp_path

        info_msgs = []
        monkeypatch.setattr("swanlab.sdk.internal.run.progress._POLL_INTERVAL", 0)
        monkeypatch.setattr(
            greeting.console, "info", lambda *args, **kwargs: info_msgs.append(" ".join(map(str, args)))
        )
        monkeypatch.setattr("swanlab.sdk.internal.run.client.exists", lambda: False)

        run.finish()

        text = "\n".join(info_msgs)
        assert "Not connected to SwanLab online service yet" in text
        assert "https://swanlab.cn/@alice/proj/runs/run-9" in text
        assert run.path == "/alice/proj/run-9"

    def test_finish_does_not_print_link_on_permanent_attach_failure(self, monkeypatch, tmp_path):
        core = MagicMock()
        core.confirm_after_deliver_failure = True
        core.last_error = None
        core.deliver_run_finish.return_value = DeliverRunFinishResponse(success=True, message="OK")
        core.confirm_run_finish.return_value = ConfirmRunFinishResponse(
            success=False, message="Failed to finish run, but it has been saved locally. 401 invalid token"
        )
        core.get_operation_stats.return_value = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(
                state=CoreState.CORE_STATE_FINISHED,
                cloud=CloudStatus(state=CLOUD_STATE_FAILED, error="401 invalid token"),
            ),
        )

        run = _make_online_run(core, path=None)
        ctx = cast(Any, run._ctx)
        ctx.run_dir = tmp_path
        ctx.config.settings = SimpleNamespace(
            mode="online",
            web_host="https://swanlab.cn",
            core=SimpleNamespace(skip_store=False),
            run=SimpleNamespace(dir=tmp_path.name),
            experiment=SimpleNamespace(name="test-exp"),
        )
        ctx.config.run_dir = tmp_path

        info_msgs = []
        warning_msgs = []
        printed_links = []
        monkeypatch.setattr(
            greeting.console, "info", lambda *args, **kwargs: info_msgs.append(" ".join(map(str, args)))
        )
        monkeypatch.setattr(
            greeting.console, "warning", lambda *args, **kwargs: warning_msgs.append(" ".join(map(str, args)))
        )
        monkeypatch.setattr(greeting, "print_cloud_link", lambda r: printed_links.append(r))
        monkeypatch.setattr("swanlab.sdk.internal.run.console.error", MagicMock())
        monkeypatch.setattr("swanlab.sdk.internal.run.client.exists", lambda: False)

        run.finish()

        assert printed_links == []
        assert "recovering in background" not in "\n".join(info_msgs)
        assert any("Failed to connect to SwanLab online service" in message for message in warning_msgs)
        assert all("swanlab sync" not in message for message in info_msgs)
        assert sum("swanlab sync" in message for message in warning_msgs) == 1
        with pytest.raises(ValueError, match="401 invalid token.*Data is saved at.*after resolving the issue"):
            _ = run.path


class TestRuntimeStatusMonitor:
    """运行期状态监听：延迟接入恢复/失败提示与运行期上传暂停/恢复提示。"""

    def test_cloud_attach_recovery_prints_link_once(self, monkeypatch):
        core = MagicMock()
        run = _make_online_run(core, path=None)

        infos: list[str] = []
        links: list[Run] = []
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.info", lambda *a, **kw: infos.append(" ".join(map(str, a)))
        )
        monkeypatch.setattr(greeting, "print_cloud_link", lambda r: links.append(r))

        pending = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(
                state=CoreState.CORE_STATE_RUNNING,
                cloud=CloudStatus(state=CLOUD_STATE_PENDING),
                upload=UploadStatus(paused=True),
            ),
        )
        attached = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(
                state=CoreState.CORE_STATE_RUNNING,
                cloud=CloudStatus(state=CLOUD_STATE_ATTACHED, path="/alice/proj/run-rec"),
                upload=UploadStatus(paused=False),
            ),
        )

        run._handle_runtime_stats(pending)
        assert links == []
        assert infos == []

        run._handle_runtime_stats(attached)
        run._handle_runtime_stats(attached)
        assert links == [run]
        assert len([m for m in infos if "Connected to SwanLab online service" in m]) == 1

    def test_cloud_attach_failed_prints_error_and_sync_once(self, monkeypatch):
        core = MagicMock()
        run = _make_online_run(core, path=None)

        errors: list[str] = []
        warnings: list[str] = []
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.error", lambda *a, **kw: errors.append(" ".join(map(str, a)))
        )
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.warning", lambda *a, **kw: warnings.append(" ".join(map(str, a)))
        )

        failed = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(
                state=CoreState.CORE_STATE_RUNNING,
                cloud=CloudStatus(state=CLOUD_STATE_FAILED, error="401 invalid token"),
                upload=UploadStatus(paused=True),
            ),
        )
        run._handle_runtime_stats(failed)
        run._handle_runtime_stats(failed)

        assert len([m for m in errors if "401 invalid token" in m]) == 1
        assert len([m for m in warnings if "swanlab sync" in m]) == 1
        assert not any("Upload paused" in m for m in warnings)

    @pytest.mark.parametrize("cloud_state", [CLOUD_STATE_UNSPECIFIED, CLOUD_STATE_ATTACHED])
    def test_runtime_upload_pause_and_resume_transitions(self, monkeypatch, cloud_state):
        core = MagicMock()
        run = _make_online_run(core, path="/alice/proj/run-1")

        warnings: list[str] = []
        infos: list[str] = []
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.warning", lambda *a, **kw: warnings.append(" ".join(map(str, a)))
        )
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.info", lambda *a, **kw: infos.append(" ".join(map(str, a)))
        )

        cloud = CloudStatus(state=cloud_state, path="/alice/proj/run-1") if cloud_state else None
        paused_resp = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(state=CoreState.CORE_STATE_RUNNING, cloud=cloud, upload=UploadStatus(paused=True)),
        )
        resumed_resp = GetOperationStatsResponse(
            success=True,
            stats=OperationStats(state=CoreState.CORE_STATE_RUNNING, cloud=cloud, upload=UploadStatus(paused=False)),
        )

        run._handle_runtime_stats(paused_resp)
        run._handle_runtime_stats(paused_resp)
        assert len([m for m in warnings if "🚨 Upload paused" in m]) == 1

        run._handle_runtime_stats(resumed_resp)
        run._handle_runtime_stats(resumed_resp)
        assert len([m for m in infos if "Upload recovered, resuming normally." in m]) == 1

    def test_monitor_thread_starts_only_for_reporting_online_core(self, monkeypatch):
        polled = threading.Event()
        core = MagicMock()
        core.reports_runtime_stats = True
        core.get_operation_stats.side_effect = lambda: (
            polled.set(),
            GetOperationStatsResponse(success=True, stats=OperationStats(state=CoreState.CORE_STATE_RUNNING)),
        )[1]

        run = _make_online_run(core, path="/alice/proj/run-1")
        monkeypatch.setattr(run, "_STATUS_POLL_INTERVAL", 0.02)
        run._start_status_monitor()
        try:
            assert run._status_timer is not None
            assert run._status_timer.is_running
            assert polled.wait(timeout=1.0)
        finally:
            run._stop_status_monitor()
        assert run._status_timer is None

        # 未显式启用 reports_runtime_stats 的 mock/Python Core 不启动监听定时器
        non_reporting_core = MagicMock()
        run_no_report = _make_online_run(non_reporting_core, path="/alice/proj/run-1")
        run_no_report._start_status_monitor()
        assert run_no_report._status_timer is None
