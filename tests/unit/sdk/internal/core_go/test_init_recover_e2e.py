"""跨进程 Python SDK <-> Go Core 断网恢复端到端测试。"""

from __future__ import annotations

import atexit
import socket
import threading
import time
from pathlib import Path
from typing import Any, Dict, Optional

import pytest

import swanlab
from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.core_python import client as py_client
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import clear_run
from swanlab.sdk.internal.settings import Settings, set_global_settings


def _base_settings(
    url: str,
    log_dir: Path,
    *,
    skip_store: bool = False,
    enable_probes: bool = False,
    extra: Optional[Dict[str, Any]] = None,
) -> Settings:
    data: Dict[str, Any] = {
        "api_key": "test-api-key",
        "api_host": url,
        "web_host": url,
        "log_dir": str(log_dir),
        "core": {
            "skip_store": skip_store,
            "health_probe_interval": 0.2,
            "health_probe_timeout": 0.5,
            "health_degrade_threshold": 1,
            "health_recover_threshold": 1,
        },
        "probe": {
            "hardware": enable_probes,
            "requirements": enable_probes,
            "conda": enable_probes,
            "git": enable_probes,
            "swanlab": enable_probes,
            "monitor": enable_probes,
        },
    }
    if extra:
        data.update(extra)
    return Settings.model_validate(data)


def test_online_init_502_recovers_and_keeps_url_after_context_exit(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """init 时后端返回 502，训练继续写入本地；退出 with 后恢复建实验并保留可读的 run.url。"""
    if py_client.exists():
        py_client.reset()
    with start_mock_cloud(server_status=502, house_status=502) as cloud:
        state = cloud.state
        with swanlab.init(
            mode="online",
            on_online_error="offline",
            project="recover-e2e",
            settings=_base_settings(cloud.url, tmp_path / "logs"),
        ) as run:
            assert run.mode == "online"
            with pytest.raises(ValueError, match="has not connected to SwanLab online service yet"):
                _ = run.path
            for step in range(1, 11):
                run.log({"loss": round(1.0 / step, 4)}, step=step)
            with state.lock:
                state.server_status = 200
                state.house_status = 200

        assert run.path == f"/@alice/recover-e2e/{run.id}"
        assert run.url == f"{cloud.url}/@alice/recover-e2e/runs/{run.id}"

    assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert sorted(state.scalar_steps) == list(range(1, 11))
    assert state.stop_states == ["FINISHED"]


def test_online_zero_record_finish_waits_for_cloud_recovery(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """未记录任何指标时调用 finish，必须等待云端恢复并完成建实验。"""
    if py_client.exists():
        py_client.reset()
    with start_mock_cloud(server_status=502, house_status=502) as cloud:
        state = cloud.state

        def _recover_later() -> None:
            time.sleep(0.6)
            with state.lock:
                state.server_status = 200
                state.house_status = 200

        timer = threading.Thread(target=_recover_later, daemon=True)
        try:
            run = swanlab.init(
                mode="online",
                on_online_error="offline",
                project="zero-wait",
                settings=_base_settings(cloud.url, tmp_path / "logs"),
            )
            t0 = time.monotonic()
            timer.start()
            run.finish()
            elapsed = time.monotonic() - t0
            assert elapsed >= 0.5
            assert run.url == f"{cloud.url}/@alice/zero-wait/runs/{run.id}"
        finally:
            timer.join(timeout=2.0)

    assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert state.stop_states == ["FINISHED"]


def test_online_skip_store_disabled_during_outage_recovery(
    tmp_path: Path, monkeypatch, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """skip_store=True 遇到恢复场景时应自动落盘并在云端恢复后补传媒体文件。"""
    if py_client.exists():
        py_client.reset()
    warnings: list[str] = []
    monkeypatch.setattr(console, "warning", lambda *args, **kwargs: warnings.append(" ".join(map(str, args))))

    artifact = tmp_path / "ckpt.txt"
    artifact.write_text("checkpoint-payload", encoding="utf-8")

    with start_mock_cloud(server_status=502, house_status=502) as cloud:
        state = cloud.state
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="skip-store-recover",
            settings=_base_settings(cloud.url, tmp_path / "logs", skip_store=True),
        )
        run.log({"acc": 0.9}, step=1)
        run.save(str(artifact), policy="now")
        assert list(run.dir.glob("run-*.swanlab")) != []
        with state.lock:
            state.server_status = 200
            state.house_status = 200
        run.finish()

    assert_no_owned_core_process()
    assert any("skip_store" in w for w in warnings)
    assert state.exp_creates == 1
    assert state.uploaded_files >= 1
    assert state.stop_states == ["FINISHED"]


def test_explicit_offline_never_contacts_healthy_cloud(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """显式 offline 模式即使后端健康也绝不访问云端。"""
    with start_mock_cloud(server_status=200, house_status=200) as cloud:
        state = cloud.state
        run = swanlab.init(
            mode="offline",
            project="strict-offline",
            settings=_base_settings(cloud.url, tmp_path / "logs"),
        )
        run.log({"loss": 0.5}, step=1)
        run.finish()

    assert_no_owned_core_process()
    assert state.total_http_requests == 0
    assert state.exp_creates == 0


def test_auth_401_degrades_offline_and_attach_three_401_reports_sync(
    tmp_path: Path, monkeypatch, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """init 直遇 401 时按策略降级 offline（raise 时抛错）；先 502 后连续 3 次 401 时提示 swanlab sync。"""
    if py_client.exists():
        py_client.reset()

    warnings: list[str] = []
    errors: list[str] = []
    monkeypatch.setattr(console, "warning", lambda *args, **kwargs: warnings.append(" ".join(map(str, args))))
    monkeypatch.setattr(console, "error", lambda *args, **kwargs: errors.append(" ".join(map(str, args))))

    # 子场景 1a：init 直接遇 401，默认降级为 offline 继续记录并在 finish 打印原因与 sync 提示。
    with start_mock_cloud(server_status=200, login_status=401) as cloud1:
        state1 = cloud1.state
        run1 = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="auth-fallback",
            settings=_base_settings(cloud1.url, tmp_path / "logs1"),
        )
        assert run1.mode == "offline"
        run1.log({"loss": 0.4}, step=1)
        run1.finish()
    assert state1.exp_creates == 0
    assert any("Authentication with SwanLab online service failed" in msg for msg in warnings)
    assert any("Authentication with SwanLab online service failed" in msg for msg in errors)
    assert_no_owned_core_process()

    # 子场景 1b：显式配置 on_online_error="raise" 时，401 立即抛 AuthenticationError。
    if py_client.exists():
        py_client.reset()
    with start_mock_cloud(server_status=200, login_status=401) as cloud1b:
        state1b = cloud1b.state
        with pytest.raises(AuthenticationError):
            swanlab.init(
                mode="online",
                on_online_error="raise",
                project="auth-raise",
                settings=_base_settings(cloud1b.url, tmp_path / "logs1b"),
            )
    assert state1b.exp_creates == 0
    assert_no_owned_core_process()

    # 子场景 2：init 时 502 进入恢复，随后健康探测 200 但认证持续 401 达 3 次，转为永久失败并提示 sync。
    if py_client.exists():
        py_client.reset()
    warnings.clear()
    errors.clear()

    with start_mock_cloud(server_status=502, house_status=502) as cloud2:
        state2 = cloud2.state
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="auth-perm-fail",
            settings=_base_settings(cloud2.url, tmp_path / "logs2"),
        )
        with state2.lock:
            state2.server_status = 200
            state2.house_status = 200
            state2.login_status = 401

        deadline = time.monotonic() + 8.0
        while time.monotonic() < deadline:
            with state2.lock:
                if state2.auth_verify_calls >= 3:
                    break
            time.sleep(0.05)
        assert state2.auth_verify_calls >= 3

        run.log({"loss": 0.3}, step=1)
        run.finish()
        with pytest.raises(ValueError, match="bad api key"):
            _ = run.path

    assert_no_owned_core_process()
    assert any("swanlab sync" in msg for msg in warnings)
    assert errors != []


def test_resume_without_local_file_falls_back_to_offline(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """resume='allow' 但目录无既有本地 run 文件时退回 offline，不进入在线恢复等待。"""
    if py_client.exists():
        py_client.reset()
    with start_mock_cloud(server_status=502, house_status=502) as cloud:
        state = cloud.state
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            resume="allow",
            id="missing-local-run",
            project="resume-no-file",
            settings=_base_settings(cloud.url, tmp_path / "logs"),
        )
        assert run.mode == "offline"
        with state.lock:
            state.server_status = 200
            state.house_status = 200
        run.log({"loss": 0.2}, step=1)
        run.finish()

    assert_no_owned_core_process()
    assert state.exp_creates == 0


def test_resume_with_existing_local_run_recovers_online_and_disables_probes(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """已有未完成本地记录且后端暂不可达时，以 resume 恢复在线并在恢复后禁用探针采集。"""
    if py_client.exists():
        py_client.reset()

    log_dir = tmp_path / "logs"
    seed_run = swanlab.init(
        mode="offline",
        project="resume-recover",
        settings=Settings.model_validate({"log_dir": str(log_dir), "probe": {"monitor": False}}),
    )
    run_id = seed_run.id
    run_dir = Path(seed_run.dir)
    for step in range(1, 4):
        swanlab.log({"loss": float(step)}, step=step)
    # 先停掉前端消费线程让离线标量真正落盘，再强行中止 Core 以模拟未 finish 的断电现场
    seed_run._components.stop()
    seed_core = seed_run._ctx.core
    if seed_run._atexit_registered:
        atexit.unregister(seed_run._handle_atexit)
        seed_run._atexit_registered = False
    clear_run()
    seed_core.abort()

    with start_mock_cloud(server_status=502, house_status=502) as cloud:
        state = cloud.state
        resumed = swanlab.init(
            mode="online",
            on_online_error="offline",
            resume="allow",
            id=run_id,
            project="resume-recover",
            settings=_base_settings(
                cloud.url,
                log_dir,
                enable_probes=True,
                extra={"run": {"dir": run_dir.name}},
            ),
        )
        assert resumed.mode == "online"
        # 既有本地未完成记录在启动恢复时 new_experiment=False，应禁用重复硬件探针
        assert resumed._ctx.config.settings.probe.hardware is False
        assert resumed._ctx.config.settings.probe.monitor is False

        with state.lock:
            state.server_status = 200
            state.house_status = 200
        resumed.log({"loss": 4.0}, step=4)
        resumed.finish()

    assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert sorted(state.scalar_steps) == [1, 2, 3, 4]
    assert state.stop_states == ["FINISHED"]


def test_explicit_login_outage_then_init_recovers_online(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    """先调用 swanlab.login() 且端口拒绝连接，不阻塞脚本；随后 init 进入后台恢复并在后端上线后完成补传。"""
    if py_client.exists():
        py_client.reset()

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe_sock:
        probe_sock.bind(("127.0.0.1", 0))
        free_port = int(probe_sock.getsockname()[1])
    target_url = f"http://127.0.0.1:{free_port}"

    cloud = None
    try:
        # 1. 目标服务未监听，显式 swanlab.login() 不抛错，返回 False 并保留目标凭证/地址
        logged_in = swanlab.login(api_key="test-api-key", host=target_url, relogin=True, timeout=1)
        assert logged_in is False

        # 2. init 在未连接状态下启动 Go Core 并进入后台恢复（直接沿用 login 写入全局配置的凭证与 host）
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="booter-recover",
            log_dir=str(tmp_path / "logs"),
            settings=Settings.model_validate(
                {
                    "probe": {
                        "hardware": False,
                        "requirements": False,
                        "conda": False,
                        "git": False,
                        "swanlab": False,
                        "monitor": False,
                    },
                    "core": {
                        "health_probe_interval": 0.2,
                        "health_probe_timeout": 0.5,
                        "health_degrade_threshold": 1,
                        "health_recover_threshold": 1,
                    },
                }
            ),
        )
        assert run.mode == "online"
        run.log({"loss": 0.5}, step=1)

        # 3. 后端恢复上线，等待 Go Core 完成附着与补传
        cloud = start_mock_cloud(server_status=200, house_status=200, port=free_port)
        state = cloud.state
        run.log({"loss": 0.2}, step=2)
        run.finish()

        assert state.exp_creates == 1
        assert sorted(state.scalar_steps) == [1, 2]
        assert state.stop_states == ["FINISHED"]
        assert run.path == f"/@alice/booter-recover/{run.id}"
        assert run.url == f"{target_url}/@alice/booter-recover/runs/{run.id}"
    finally:
        if cloud is not None:
            cloud.close()
        if py_client.exists():
            py_client.reset()
        set_global_settings(Settings())

    assert_no_owned_core_process()
