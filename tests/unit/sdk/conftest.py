"""SDK 跨进程端到端测试共享 Mock 云端与 Core 进程校验辅助。"""

from __future__ import annotations

import json
import os
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import pytest

from swanlab.sdk.internal.core_go.config import CoreGoError
from swanlab.sdk.internal.core_go.runtime import locate_core_binary


class MockCloudState:
    """记录模拟 SwanLab 云端收到的请求状态，并支持动态切换响应状态码。"""

    def __init__(
        self,
        *,
        server_status: int = 200,
        house_status: int = 200,
        login_status: Optional[int] = None,
        existing_cuids: Optional[Dict[str, str]] = None,
        summary: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.lock = threading.Lock()
        self.server_status = server_status
        self.house_status = house_status
        self.login_status = login_status
        self.summary: Dict[str, Any] = dict(summary) if summary is not None else {"scalar": []}
        self.total_http_requests = 0
        self.auth_verify_calls = 0
        self.exp_creates = 0
        self.exp_requests: List[Dict[str, Any]] = []
        self.created_cuids: Dict[str, str] = dict(existing_cuids or {})
        self.last_create_body: Dict[str, Any] = {}
        self.series_keys: List[str] = []
        self.stop_states: List[str] = []
        self.profiles: List[Dict[str, Any]] = []
        self.scalar_steps: List[int] = []
        self.scalar_steps_by_exp: Dict[str, List[int]] = {}
        self.prepared_files: List[str] = []
        self.uploaded_files = 0
        self.completed_files: List[str] = []


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True


class MockCloudServer:
    """封装本地线程化 HTTP Mock 服务及其幂等关闭生命周期。"""

    def __init__(self, server: ThreadingHTTPServer, state: MockCloudState, url: str) -> None:
        self.server = server
        self.state = state
        self.url = url
        self._closed = False

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self.server.shutdown()
        self.server.server_close()

    def __enter__(self) -> MockCloudServer:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def _create_mock_cloud(
    *,
    server_status: int = 200,
    house_status: int = 200,
    login_status: Optional[int] = None,
    port: int = 0,
    existing_cuids: Optional[Dict[str, str]] = None,
    summary: Optional[Dict[str, Any]] = None,
    use_request_cuid: bool = True,
) -> MockCloudServer:
    state = MockCloudState(
        server_status=server_status,
        house_status=house_status,
        login_status=login_status,
        existing_cuids=existing_cuids,
        summary=summary,
    )
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
                    state.exp_requests.append(body)
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
                    exp_id = cuid if (use_request_cuid and cuid) else f"exp-{state.exp_creates}"
                    if use_request_cuid and cuid:
                        state.created_cuids[cuid] = exp_id
                self._write_json(201, {"cuid": exp_id, "slug": exp_id, "name": body.get("name", "")})
                return

            if len(parts) == 5 and parts[:2] == ["api", "projects"] and parts[4] == "series":
                code = self._server_code()
                body = self._read_json()
                if code == 200:
                    keys = [item.get("key") for item in body.get("series") or [] if isinstance(item, dict)]
                    with state.lock:
                        state.series_keys.extend(str(k) for k in keys if k is not None)
                self._write_json(code, {"ok": code == 200})
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
                exp_id = str(body.get("experimentId") or "")
                if body.get("type") == "scalar":
                    with state.lock:
                        steps_for_exp = state.scalar_steps_by_exp.setdefault(exp_id, [])
                        for item in body.get("metrics") or []:
                            if isinstance(item, dict) and isinstance(item.get("index"), int):
                                idx = item["index"]
                                state.scalar_steps.append(idx)
                                steps_for_exp.append(idx)
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
                self._write_json(200, state.summary)
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

            if len(parts) == 7 and parts[:2] == ["api", "project"] and parts[4] == "runs" and parts[6] == "profile":
                body = self._read_json()
                with state.lock:
                    state.profiles.append(body)
                self._write_json(200, {"ok": True})
                return

            self._read_json()
            self._write_json(200, {"ok": True})

    server = ReusableThreadingHTTPServer(("127.0.0.1", port), Handler)
    host, bound_port = server.server_address[:2]
    url = f"http://{host!s}:{bound_port}"
    server_ref["url"] = url
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return MockCloudServer(server, state, url)


@pytest.fixture
def require_core_binary() -> Path:
    """确认本机存在 swanlab-core 可执行文件；缺失时跳过端到端测试。"""
    try:
        return locate_core_binary()
    except CoreGoError as exc:
        pytest.skip(f"swanlab-core binary unavailable: {exc}")


@pytest.fixture
def assert_no_owned_core_process() -> Callable[[], None]:
    """返回检查当前测试进程名下无残留 swanlab-core 子进程的断言函数。"""

    def _assert() -> None:
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

    return _assert


@pytest.fixture
def start_mock_cloud():
    """启动本地模拟云端服务器并在测试结束时自动关闭所有实例。"""
    servers: List[MockCloudServer] = []

    def _start(
        *,
        server_status: int = 200,
        house_status: int = 200,
        login_status: Optional[int] = None,
        port: int = 0,
        existing_cuids: Optional[Dict[str, str]] = None,
        summary: Optional[Dict[str, Any]] = None,
        use_request_cuid: bool = True,
    ) -> MockCloudServer:
        srv = _create_mock_cloud(
            server_status=server_status,
            house_status=house_status,
            login_status=login_status,
            port=port,
            existing_cuids=existing_cuids,
            summary=summary,
            use_request_cuid=use_request_cuid,
        )
        servers.append(srv)
        return srv

    try:
        yield _start
    finally:
        for srv in reversed(servers):
            srv.close()
