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
    CloudStatus,
    CoreState,
    OperationStats,
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

        with pytest.raises(ValueError, match="has not attached yet"):
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
        assert "SwanLab Cloud has not attached yet" in text
        assert "https://swanlab.cn/@alice/proj/runs/run-9" in text
        assert run.path == "/alice/proj/run-9"

    def test_finish_does_not_print_link_on_permanent_attach_failure(self, monkeypatch, tmp_path):
        core = MagicMock()
        core.confirm_after_deliver_failure = True
        core.last_error = None
        core.deliver_run_finish.return_value = DeliverRunFinishResponse(success=True, message="OK")
        core.confirm_run_finish.return_value = ConfirmRunFinishResponse(success=False, message="401 invalid token")
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
        assert any("Cloud recovery failed" in message for message in warning_msgs)
        assert any("swanlab sync" in message for message in info_msgs)
        with pytest.raises(ValueError, match="401 invalid token.*swanlab sync"):
            _ = run.path
