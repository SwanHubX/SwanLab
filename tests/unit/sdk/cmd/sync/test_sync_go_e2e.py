"""
@author: Nexisato
@file: test_sync_go_e2e.py
@time: 2026/11/18
@description: 默认 Go core 路由下 swanlab.init(offline) -> swanlab.sync() 端到端集成测试。
"""

from __future__ import annotations

import json
import os
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pytest

import swanlab
from swanlab.sdk.cmd.login import login_raw
from swanlab.sdk.cmd.sync import sync as sync_run
from swanlab.sdk.internal import impl
from swanlab.sdk.internal.core_go.config import CoreGoError
from swanlab.sdk.internal.core_go.runtime import locate_core_binary
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings.experiment import ProjectSettings, RunSettings


class _MockCloudState:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.experiments = 0
        self.exp_requests: List[Dict[str, Any]] = []
        self.stop_states: List[str] = []
        self.scalar_steps: Dict[str, List[int]] = {}
        self.series_keys: List[str] = []
        self.profiles: List[Dict[str, Any]] = []


def _start_mock_cloud() -> Tuple[ThreadingHTTPServer, _MockCloudState, str]:
    state = _MockCloudState()

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

        def _read_json(self) -> Dict[str, Any]:
            length = int(self.headers.get("Content-Length", "0") or 0)
            raw = self.rfile.read(length) if length > 0 else b"{}"
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

        def do_POST(self) -> None:
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]
            if parts == ["api", "login", "api_key"]:
                self._read_json()
                self._write_json(
                    200,
                    {
                        "sid": "mock-sid",
                        "expiredAt": "2099-01-01T00:00:00.000Z",
                        "userInfo": {"username": "alice"},
                    },
                )
                return
            if len(parts) == 5 and parts[:2] == ["api", "project"] and parts[4] == "experiment":
                body = self._read_json()
                with state.lock:
                    state.exp_requests.append(body)
                if body.get("cuid") == "existing-run":
                    self._write_json(
                        200,
                        {
                            "cuid": "existing-run",
                            "slug": "existing-run",
                            "name": body.get("name", ""),
                            "createdAt": "2026-01-02T00:00:00.000Z",
                        },
                    )
                    return
                with state.lock:
                    state.experiments += 1
                    exp_id = f"exp-{state.experiments}"
                self._write_json(201, {"cuid": exp_id, "slug": exp_id, "name": body.get("name", "")})
                return
            if len(parts) == 5 and parts[:2] == ["api", "projects"] and parts[4] == "series":
                body = self._read_json()
                keys = [item.get("key") for item in body.get("series") or [] if isinstance(item, dict)]
                with state.lock:
                    state.series_keys.extend(str(k) for k in keys if k is not None)
                self._write_json(200, {"ok": True})
                return
            if parts == ["api", "house", "metrics"]:
                body = self._read_json()
                exp_id = str(body.get("experimentId") or "")
                metric_type = str(body.get("type") or "")
                items = body.get("metrics") or []
                if metric_type == "scalar":
                    with state.lock:
                        steps = state.scalar_steps.setdefault(exp_id, [])
                        for item in items:
                            if isinstance(item, dict) and isinstance(item.get("index"), int):
                                steps.append(item["index"])
                self._write_json(200, {"ok": True})
                return
            self._write_json(200, {"ok": True})

        def do_GET(self) -> None:
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]
            if parts == ["api", "auth", "verify"]:
                self._write_json(200, {"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
                return
            if len(parts) == 4 and parts[:2] == ["api", "project"]:
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
            if parts[:4] == ["api", "house", "metrics", "summaries"]:
                self._write_json(200, {"scalar": [{"key": "loss", "step": 12}]})
                return
            self._write_json(404, {"code": "not_found"})

        def do_PUT(self) -> None:
            path = self.path.split("?", 1)[0].rstrip("/")
            parts = [p for p in path.split("/") if p]
            if len(parts) == 7 and parts[:2] == ["api", "project"] and parts[4] == "runs" and parts[6] == "state":
                body = self._read_json()
                with state.lock:
                    if isinstance(body.get("state"), str):
                        state.stop_states.append(body["state"])
                self._write_json(200, {"ok": True})
                return
            if len(parts) == 7 and parts[:2] == ["api", "project"] and parts[4] == "runs" and parts[6] == "profile":
                body = self._read_json()
                with state.lock:
                    state.profiles.append(body)
                self._write_json(200, {"ok": True})
                return
            self._write_json(200, {"ok": True})

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address[:2]
    return server, state, f"http://{host!s}:{port}"


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


def _create_offline_run_dir(*, steps: int = 5, finish_run: bool = True, log_dir: Optional[Path] = None) -> Path:
    settings = Settings(log_dir=log_dir) if log_dir is not None else Settings()
    run = swanlab.init(mode="offline", project="sync-go-e2e", settings=settings)
    run_dir = Path(run._ctx.run_dir)
    for step in range(1, steps + 1):
        swanlab.log({"loss": round(1.0 / step, 4)}, step=step)
    if finish_run:
        run.finish()
    else:
        from swanlab.sdk.internal.run import clear_run

        core = run._ctx.core
        clear_run()
        core.abort()
    return run_dir


def test_sync_default_go_full_with_active_client_api_key_fallback(tmp_path: Path):
    _require_core_binary()
    assert impl.resolve_core_sync_backend().impl == "go"

    run_dir = _create_offline_run_dir(steps=5, log_dir=tmp_path / "logs")
    server, state, url = _start_mock_cloud()
    try:
        login_raw("session-only-key", host=url, print_welcome=False)
        sync_run(
            run_dir,
            settings=Settings(api_host=url, web_host=url, project=ProjectSettings(name="sync-go-e2e")),
        )
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
    assert state.experiments == 1
    assert state.stop_states == ["FINISHED"]
    assert sorted(state.scalar_steps["exp-1"]) == [1, 2, 3, 4, 5]
    assert "loss" in state.series_keys


def test_sync_default_go_incremental_range_and_corrupt_tail(tmp_path: Path):
    _require_core_binary()
    run_dir = _create_offline_run_dir(steps=20, log_dir=tmp_path / "logs")
    server, state, url = _start_mock_cloud()
    try:
        sync_settings = Settings(
            api_host=url,
            web_host=url,
            api_key="smoke-key",
            project=ProjectSettings(name="sync-go-e2e"),
            run=RunSettings(id="existing-run"),
        )
        # 1. 增量同步：mock summary 为 step=12，仅应补传 13..20
        sync_run(run_dir, settings=sync_settings)
        assert sorted(state.scalar_steps["existing-run"]) == list(range(13, 21))

        # 2. 范围回放：start_step=5, end_step=10，忽略 summary 水位 12 重放 5..10
        state.scalar_steps.clear()
        sync_run(run_dir, settings=sync_settings, start_step=5, end_step=10)
        assert sorted(state.scalar_steps["existing-run"]) == list(range(5, 11))

        # 3. 损坏尾部容错：截断末尾 4 字节后仍可同步前缀记录
        run_files = list(run_dir.glob("run-*.swanlab"))
        assert len(run_files) == 1
        raw = run_files[0].read_bytes()
        run_files[0].write_bytes(raw[:-4])

        corrupt_dir = tmp_path / "corrupt-copy"
        corrupt_dir.mkdir(parents=True, exist_ok=True)
        for p in run_dir.iterdir():
            if p.is_file():
                (corrupt_dir / p.name).write_bytes(p.read_bytes())

        state.scalar_steps.clear()
        state.stop_states.clear()
        sync_run(
            corrupt_dir,
            settings=Settings(
                api_host=url, web_host=url, api_key="smoke-key", project=ProjectSettings(name="sync-go-e2e")
            ),
        )
        uploaded = state.scalar_steps["exp-1"]
        assert len(uploaded) > 0
        assert uploaded[0] == 1
        assert len(uploaded) <= 20
        # 截断掉了尾部 FinishRecord，新实验按 CRASHED 收尾
        assert state.stop_states == ["CRASHED"]
    finally:
        server.shutdown()
        server.server_close()

    _assert_no_owned_core_process()
