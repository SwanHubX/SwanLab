"""跨进程 Python SDK <-> Go Core 断网恢复端到端测试。"""

from __future__ import annotations

import json
import os
import subprocess
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pytest

import swanlab
from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.core_go.config import CoreGoError
from swanlab.sdk.internal.core_go.runtime import locate_core_binary
from swanlab.sdk.internal.core_python import client as py_client
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import clear_run
from swanlab.sdk.internal.settings import Settings


class _FaultyCloudState:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.server_status = 200
        self.house_status = 200
        self.login_status: Optional[int] = None
        self.total_http_requests = 0
        self.auth_verify_calls = 0
        self.exp_creates = 0
        self.created_cuids: Dict[str, str] = {}
        self.stop_states: List[str] = []
        self.scalar_steps: List[int] = []
        self.prepared_files: List[str] = []
        self.uploaded_files = 0
        self.completed_files: List[str] = []
        self.last_create_body: Dict[str, Any] = {}


def _start_faulty_cloud(
    *,
    server_status: int = 200,
    house_status: int = 200,
    login_status: Optional[int] = None,
) -> Tuple[ThreadingHTTPServer, _FaultyCloudState, str]:
    state = _FaultyCloudState()
    state.server_status = server_status
    state.house_status = house_status
    state.login_status = login_status

    server_ref: Dict[str, str] = {}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

        def _read_body(self) -> bytes:
            length = int(self.headers.get("Content-Length", "0") or 0)
            return self.rfile.read(length) if length > 0 else b""

        def _read_json(self) -> Dict[str, Any]:
            raw = self._read_body()
            if not raw:
                return {}
            try:
                parsed = json.loads(raw)
                return parsed if isinstance(parsed, dict) else {}
            except json.JSONDecodeError:
                return {}

        def _write_json(self, status: int, payload: Dict[str, Any]) -> None:
            data = json.dumps(payload).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _server_code(self, *, auth: bool = False) -> int:
            with state.lock:
                if auth and state.login_status is not None:
                    return state.login_status
                return state.server_status

        def _house_code(self) -> int:
            with state.lock:
                return state.house_status

        def do_POST(self) -> None:
            with state.lock:
                state.total_http_requests += 1
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]

            if parts == ["api", "login", "api_key"]:
                self._read_json()
                code = self._server_code(auth=True)
                if code != 200:
                    self._write_json(code, {"code": "auth_or_outage", "message": f"login status {code}"})
                    return
                self._write_json(
                    200,
                    {
                        "sid": "mock-sid",
                        "expiredAt": "2099-01-01T00:00:00.000Z",
                        "userInfo": {"username": "alice", "name": "Alice"},
                    },
                )
                return

            if len(parts) == 5 and parts[:2] == ["api", "project"] and parts[4] == "experiment":
                code = self._server_code()
                if code != 200:
                    self._write_json(code, {"message": f"experiment unavailable ({code})"})
                    return
                body = self._read_json()
                cuid = str(body.get("cuid") or "")
                with state.lock:
                    state.last_create_body = body
                    if cuid and cuid in state.created_cuids:
                        exp_id = state.created_cuids[cuid]
                        self._write_json(
                            200,
                            {
                                "cuid": exp_id,
                                "slug": exp_id,
                                "name": body.get("name", ""),
                                "createdAt": "2026-01-02T00:00:00.000Z",
                            },
                        )
                        return
                    state.exp_creates += 1
                    exp_id = cuid or f"exp-{state.exp_creates}"
                    if cuid:
                        state.created_cuids[cuid] = exp_id
                self._write_json(201, {"cuid": exp_id, "slug": exp_id, "name": body.get("name", "")})
                return

            if len(parts) == 5 and parts[:2] == ["api", "projects"] and parts[4] == "series":
                code = self._server_code()
                self._read_json()
                self._write_json(code if code != 200 else 200, {"ok": code == 200})
                return

            if len(parts) == 5 and parts[:2] == ["api", "experiment"] and parts[3:] == ["files", "prepare"]:
                body = self._read_json()
                files = body.get("files") or []
                urls: List[str] = []
                with state.lock:
                    for item in files:
                        rel = str(item.get("path") or "unknown")
                        state.prepared_files.append(rel)
                        urls.append(f"{server_ref['url']}/mock-oss/{rel}")
                self._write_json(200, {"urls": urls})
                return

            if len(parts) == 5 and parts[:2] == ["api", "experiment"] and parts[3:] == ["files", "complete"]:
                body = self._read_json()
                files = body.get("files") or []
                with state.lock:
                    for item in files:
                        if isinstance(item, dict) and item.get("path"):
                            state.completed_files.append(str(item["path"]))
                self._write_json(201, {"ok": True})
                return

            if parts == ["api", "house", "metrics"]:
                code = self._house_code()
                if code != 200:
                    self._read_json()
                    self._write_json(code, {"message": f"house status {code}"})
                    return
                body = self._read_json()
                if body.get("type") == "scalar":
                    with state.lock:
                        for item in body.get("metrics") or []:
                            if isinstance(item, dict) and isinstance(item.get("index"), int):
                                state.scalar_steps.append(item["index"])
                self._write_json(200, {"ok": True})
                return

            self._read_json()
            self._write_json(200, {"ok": True})

        def do_GET(self) -> None:
            with state.lock:
                state.total_http_requests += 1
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]

            if parts == ["api", "health"]:
                code = self._server_code()
                self._write_json(code, {"ok": code == 200})
                return

            if parts == ["api", "house", "health"]:
                code = self._house_code()
                self._write_json(code, {"ok": code == 200})
                return

            if parts == ["api", "auth", "verify"]:
                with state.lock:
                    state.auth_verify_calls += 1
                code = self._server_code(auth=True)
                if code != 200:
                    self._write_json(code, {"code": "auth_error", "message": f"bad api key ({code})"})
                    return
                self._write_json(200, {"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
                return

            if len(parts) == 4 and parts[:2] == ["api", "project"]:
                code = self._server_code()
                if code != 200:
                    self._write_json(code, {"message": "project unavailable"})
                    return
                proj = parts[3]
                self._write_json(
                    200,
                    {
                        "cuid": "proj-1",
                        "name": proj,
                        "path": f"/@alice/{proj}",
                        "visibility": "PRIVATE",
                        "version": 1,
                        "_count": {"experiments": 1},
                        "group": {"username": "alice"},
                    },
                )
                return

            if len(parts) == 6 and parts[:2] == ["api", "project"] and parts[4] == "runs":
                cuid = parts[5]
                with state.lock:
                    exists = cuid in state.created_cuids
                if exists:
                    self._write_json(200, {"cuid": cuid, "slug": cuid, "name": cuid})
                else:
                    self._write_json(404, {"code": "not_found"})
                return

            if parts[:4] == ["api", "house", "metrics", "summaries"]:
                self._write_json(200, {"scalar": []})
                return

            self._write_json(404, {"code": "not_found"})

        def do_PUT(self) -> None:
            with state.lock:
                state.total_http_requests += 1
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]

            if len(parts) >= 2 and parts[0] == "mock-oss":
                self._read_body()
                with state.lock:
                    state.uploaded_files += 1
                self.send_response(200)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return

            if len(parts) == 7 and parts[:2] == ["api", "project"] and parts[4] == "runs" and parts[6] == "state":
                body = self._read_json()
                with state.lock:
                    if isinstance(body.get("state"), str):
                        state.stop_states.append(body["state"])
                self._write_json(200, {"ok": True})
                return

            self._read_json()
            self._write_json(200, {"ok": True})

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    host, port = server.server_address[:2]
    url = f"http://{host!s}:{port}"
    server_ref["url"] = url
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, state, url


def _require_core_binary() -> Path:
    try:
        return locate_core_binary()
    except CoreGoError as exc:
        pytest.skip(f"swanlab-core binary unavailable: {exc}")


def _assert_no_owned_core_process() -> None:
    if os.name == "nt":
        return
    output = subprocess.run(
        ["ps", "-axo", "pid=,command="],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    marker = f"--parent-pid {os.getpid()}"
    leftover = [line.strip() for line in output.splitlines() if "swanlab-core" in line and marker in line]
    assert leftover == [], f"orphan swanlab-core process found: {leftover}"


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


def test_online_init_502_recovers_and_keeps_url_after_context_exit(tmp_path: Path):
    """init 时后端返回 502，训练继续写入本地；退出 with 后恢复建实验并保留可读的 run.url。"""
    _require_core_binary()
    if py_client.exists():
        py_client.reset()
    server, state, url = _start_faulty_cloud(server_status=502, house_status=502)
    try:
        with swanlab.init(
            mode="online",
            on_online_error="offline",
            project="recover-e2e",
            settings=_base_settings(url, tmp_path / "logs"),
        ) as run:
            assert run.mode == "online"
            with pytest.raises(ValueError, match="has not attached yet"):
                _ = run.path
            for step in range(1, 11):
                run.log({"loss": round(1.0 / step, 4)}, step=step)
            with state.lock:
                state.server_status = 200
                state.house_status = 200

        assert run.path == f"/@alice/recover-e2e/{run.id}"
        assert run.url == f"{url}/@alice/recover-e2e/runs/{run.id}"
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert sorted(state.scalar_steps) == list(range(1, 11))
    assert state.stop_states == ["FINISHED"]


def test_online_zero_record_finish_waits_for_cloud_recovery(tmp_path: Path):
    """未记录任何指标时调用 finish，必须等待云端恢复并完成建实验。"""
    _require_core_binary()
    if py_client.exists():
        py_client.reset()
    server, state, url = _start_faulty_cloud(server_status=502, house_status=502)

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
            settings=_base_settings(url, tmp_path / "logs"),
        )
        t0 = time.monotonic()
        timer.start()
        run.finish()
        elapsed = time.monotonic() - t0
        assert elapsed >= 0.5
        assert run.url == f"{url}/@alice/zero-wait/runs/{run.id}"
    finally:
        timer.join(timeout=2.0)
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert state.stop_states == ["FINISHED"]


def test_online_skip_store_disabled_during_outage_recovery(tmp_path: Path, monkeypatch):
    """skip_store=True 遇到恢复场景时应自动落盘并在云端恢复后补传媒体文件。"""
    _require_core_binary()
    if py_client.exists():
        py_client.reset()
    warnings: List[str] = []
    monkeypatch.setattr(console, "warning", lambda *args, **kwargs: warnings.append(" ".join(map(str, args))))

    artifact = tmp_path / "ckpt.txt"
    artifact.write_text("checkpoint-payload", encoding="utf-8")

    server, state, url = _start_faulty_cloud(server_status=502, house_status=502)
    try:
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="skip-store-recover",
            settings=_base_settings(url, tmp_path / "logs", skip_store=True),
        )
        run.log({"acc": 0.9}, step=1)
        run.save(str(artifact), policy="now")
        assert list(run.dir.glob("run-*.swanlab")) != []
        with state.lock:
            state.server_status = 200
            state.house_status = 200
        run.finish()
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert any("skip_store" in w for w in warnings)
    assert state.exp_creates == 1
    assert state.uploaded_files >= 1
    assert state.stop_states == ["FINISHED"]


def test_explicit_offline_never_contacts_healthy_cloud(tmp_path: Path):
    """显式 offline 模式即使后端健康也绝不访问云端。"""
    _require_core_binary()
    server, state, url = _start_faulty_cloud(server_status=200, house_status=200)
    try:
        run = swanlab.init(
            mode="offline",
            project="strict-offline",
            settings=_base_settings(url, tmp_path / "logs"),
        )
        run.log({"loss": 0.5}, step=1)
        run.finish()
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.total_http_requests == 0
    assert state.exp_creates == 0


def test_auth_401_fails_fast_and_attach_three_401_reports_sync(tmp_path: Path, monkeypatch):
    """401 启动直接报错；先 502 后连续 3 次 401 时本地继续运行、finish 不阻塞并提示 swanlab sync。"""
    _require_core_binary()
    if py_client.exists():
        py_client.reset()

    # 子场景 1：init 直接遇 401，立即抛 AuthenticationError。
    server1, state1, url1 = _start_faulty_cloud(server_status=200, login_status=401)
    try:
        with pytest.raises(AuthenticationError):
            swanlab.init(
                mode="online",
                on_online_error="offline",
                project="auth-fast-fail",
                settings=_base_settings(url1, tmp_path / "logs1"),
            )
    finally:
        server1.shutdown()
        server1.server_close()
    assert state1.exp_creates == 0
    _assert_no_owned_core_process()

    # 子场景 2：init 时 502 进入恢复，随后健康探测 200 但认证持续 401 达 3 次，转为永久失败并提示 sync。
    if py_client.exists():
        py_client.reset()
    warnings: List[str] = []
    errors: List[str] = []
    monkeypatch.setattr(console, "warning", lambda *args, **kwargs: warnings.append(" ".join(map(str, args))))
    monkeypatch.setattr(console, "error", lambda *args, **kwargs: errors.append(" ".join(map(str, args))))

    server2, state2, url2 = _start_faulty_cloud(server_status=502, house_status=502)
    try:
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            project="auth-perm-fail",
            settings=_base_settings(url2, tmp_path / "logs2"),
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
    finally:
        server2.shutdown()
        server2.server_close()

    _assert_no_owned_core_process()
    assert any("swanlab sync" in msg for msg in warnings)
    assert errors != []


def test_resume_without_local_file_falls_back_to_offline(tmp_path: Path):
    """resume='allow' 但目录无既有本地 run 文件时退回 offline，不进入在线恢复等待。"""
    _require_core_binary()
    if py_client.exists():
        py_client.reset()
    server, state, url = _start_faulty_cloud(server_status=502, house_status=502)
    try:
        run = swanlab.init(
            mode="online",
            on_online_error="offline",
            resume="allow",
            id="missing-local-run",
            project="resume-no-file",
            settings=_base_settings(url, tmp_path / "logs"),
        )
        assert run.mode == "offline"
        with state.lock:
            state.server_status = 200
            state.house_status = 200
        run.log({"loss": 0.2}, step=1)
        run.finish()
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.exp_creates == 0


def test_resume_with_existing_local_run_recovers_online_and_disables_probes(tmp_path: Path):
    """已有未完成本地记录且后端暂不可达时，以 resume 恢复在线并在恢复后禁用探针采集。"""
    _require_core_binary()
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
    clear_run()
    seed_core.abort()

    server, state, url = _start_faulty_cloud(server_status=502, house_status=502)
    try:
        resumed = swanlab.init(
            mode="online",
            on_online_error="offline",
            resume="allow",
            id=run_id,
            project="resume-recover",
            settings=_base_settings(
                url,
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
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert sorted(state.scalar_steps) == [1, 2, 3, 4]
    assert state.stop_states == ["FINISHED"]
