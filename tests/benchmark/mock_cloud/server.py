#!/usr/bin/env python3
"""
@author: nexisato
@file: server.py
@description: 压测用的本地 Mock Cloud：只累加计数，不保存请求体与指标明细。

移植自 playground/core-plan/io-plan/smoke/mockcloud.py，路由与响应体保持一致，
只做三处改动：

1. 不缓存明细，只按 record 类型累加条数、批次数和字节数，90 万条记录也只占几个整数；
2. 补齐 /api/health 与 /api/house/health，缺失会让 Go core 认为站点不可达并暂停上传；
3. 未注册路径返回 200 空体，5xx 或断连同样会被判为不可达。

同一个 run id 重复创建实验时返回同一个 experimentId，供多进程共享一个实验。

用法::

    python3 server.py --port 0 --port-file /tmp/mock.port [--latency-ms 20]
    curl -s http://127.0.0.1:<port>/__mock/stats | jq
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, List, Optional, Tuple, cast
from urllib.parse import urlparse

# 指标类型白名单，与 core/internal/api/model.go 的 MetricKind 一致。
METRIC_TYPES = ("scalar", "log", "media", "histogram")

# 单批条数的环形缓冲长度（只用于观察批量分布）。
BATCH_RING = 64

VERIFY_PROFILE = {
    "uid": 10001,
    "username": "bench",
    "name": "Bench Tester",
    "avatar": "",
    "verified": False,
    "createdAt": "2026-01-01T00:00:00.000Z",
}


class CloudState:
    """进程内计数状态；所有写操作持锁，读快照时拷贝。"""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        # 指标
        self.metrics_by_type: Dict[str, int] = {name: 0 for name in METRIC_TYPES}
        self.metrics_batches = 0
        self.metrics_bytes = 0
        self.batch_ring: List[int] = []
        self.first_metrics_at: Optional[float] = None
        self.last_metrics_at: Optional[float] = None
        # 列定义
        self.series_keys: List[str] = []
        self._series_seen = set()
        self.series_requests = 0
        # 项目 / 实验
        self.projects = set()
        self.experiments: Dict[str, str] = {}
        self.experiment_seq = 0
        # 心跳与终态
        self.heartbeats = 0
        self.stop_experiments = 0
        self.profile_uploads = 0
        # 文件
        self.presigned_requests = 0
        self.presigned_paths = 0
        self.files_prepared = 0
        self.files_completed = 0
        self.sink_files = 0
        self.sink_bytes = 0
        # 路由命中计数（未知路由记在 __other__）
        self.routes: Dict[str, int] = {}
        # 拒绝（非 200 的指标批次）
        self.rejected = 0

    # -- 写入 ----------------------------------------------------------------

    def count_route(self, method: str, route: str) -> None:
        key = f"{method} {route}"
        with self.lock:
            self.routes[key] = self.routes.get(key, 0) + 1

    def ingest_metrics(self, metric_type: str, metrics: List[Any], size: int, status: int) -> None:
        now = time.monotonic()
        with self.lock:
            self.metrics_batches += 1
            self.metrics_bytes += size
            self.batch_ring.append(len(metrics))
            if len(self.batch_ring) > BATCH_RING:
                del self.batch_ring[: len(self.batch_ring) - BATCH_RING]
            if self.first_metrics_at is None:
                self.first_metrics_at = now
            self.last_metrics_at = now
            if status != 200:
                self.rejected += 1
                return
            bucket = metric_type if metric_type in self.metrics_by_type else "log"
            self.metrics_by_type[bucket] += len(metrics)

    def ingest_columns(self, keys: List[str]) -> None:
        with self.lock:
            self.series_requests += 1
            for key in keys:
                if key not in self._series_seen:
                    self._series_seen.add(key)
                    self.series_keys.append(key)

    def ensure_experiment(self, run_id: str) -> Tuple[str, bool]:
        """按 run cuid 幂等分配 experimentId：返回 (id, 是否新建)。"""
        with self.lock:
            existing = self.experiments.get(run_id)
            if existing is not None:
                return existing, False
            self.experiment_seq += 1
            experiment_id = f"bench-exp-{self.experiment_seq}"
            self.experiments[run_id] = experiment_id
            return experiment_id, True

    def record_sink(self, size: int, digest: str) -> None:
        with self.lock:
            self.sink_files += 1
            self.sink_bytes += size

    # -- 读取 ----------------------------------------------------------------

    def snapshot(self) -> Dict[str, Any]:
        with self.lock:
            house_total = sum(self.metrics_by_type.values())
            return {
                "metrics_by_type": dict(self.metrics_by_type),
                "house_metrics_total": house_total,
                "metrics_batches": self.metrics_batches,
                "metrics_bytes": self.metrics_bytes,
                "batch_ring": list(self.batch_ring),
                "rejected_metrics": self.rejected,
                "first_metrics_at": self.first_metrics_at,
                "last_metrics_at": self.last_metrics_at,
                "series_keys": list(self.series_keys),
                "series_requests": self.series_requests,
                "projects": len(self.projects),
                "experiments": dict(self.experiments),
                "heartbeats": self.heartbeats,
                "stop_experiments": self.stop_experiments,
                "profile_uploads": self.profile_uploads,
                "presigned_requests": self.presigned_requests,
                "presigned_paths": self.presigned_paths,
                "files_prepared": self.files_prepared,
                "files_completed": self.files_completed,
                "sink_files": self.sink_files,
                "sink_bytes": self.sink_bytes,
                "routes": dict(self.routes),
            }


class MockHandler(BaseHTTPRequestHandler):
    """只覆盖 bench 需要的云端接口，未命中一律 200。"""

    protocol_version = "HTTP/1.1"
    state: CloudState
    latency_ms: int = 0

    # ------------------------------------------------------------ 基础设施

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002 - 与基类签名对齐
        if os.environ.get("BENCH_MOCK_VERBOSE"):
            sys.stderr.write("[mock] " + format % args + "\n")

    def _port(self) -> int:
        """当前监听端口。"""
        return int(cast(Any, self.server).server_address[1])

    def _read_body(self) -> bytes:
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = 0
        if length <= 0:
            return b""
        return self.rfile.read(length)

    def _json(self, status: int, payload: Any) -> None:
        if self.latency_ms:
            time.sleep(self.latency_ms / 1000.0)
        body = json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _text(self, status: int, text: str) -> None:
        body = text.encode()
        self.send_response(status)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802 - http.server 约定
        self._dispatch("GET")

    def do_POST(self) -> None:  # noqa: N802
        self._dispatch("POST")

    def do_PUT(self) -> None:  # noqa: N802
        self._dispatch("PUT")

    # ------------------------------------------------------------ 路由

    def _dispatch(self, method: str) -> None:
        parsed = urlparse(self.path)
        path = parsed.path
        body = self._read_body()

        if path.startswith("/__mock/"):
            self._control(path, body)
            return
        if path.startswith("/sink/"):
            self._sink(path, body)
            return
        if not path.startswith("/api/"):
            # 非 /api 路径直接返回 200。
            self._json(200, {"ok": True})
            return

        route = path[len("/api") :]
        self.state.count_route(method, route)
        handler = self._match(method, route)
        if handler is None:
            # 未注册路径统一返回 200，5xx 或断连会让 core 认为站点不可达。
            self._json(200, {"ok": True})
            return
        handler(body, route)

    def _control(self, path: str, body: bytes) -> None:
        if path == "/__mock/stats":
            self._json(200, self.state.snapshot())
            return
        if path == "/__mock/health":
            self._json(200, {"ok": True})
            return
        if path == "/__mock/shutdown":
            self._json(200, {"ok": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        self._json(404, {"code": "not_found", "message": path})

    def _sink(self, path: str, body: bytes) -> None:
        """接收对象存储上传，只统计数量与体积，丢弃内容。"""
        self.state.record_sink(len(body), hashlib.md5(body).hexdigest())
        self.send_response(200)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _match(self, method: str, route: str):
        parts = [p for p in route.split("/") if p]

        # 认证与健康
        if method == "GET" and route == "/auth/verify":
            return self._verify
        if method == "GET" and route == "/health":
            return self._health
        if method == "GET" and route == "/house/health":
            return self._health
        if method == "POST" and route == "/login/api_key":
            return self._login

        # 项目
        if method == "GET" and len(parts) == 3 and parts[0] == "project":
            return self._get_project
        if method == "POST" and len(parts) == 2 and parts[0] == "projects":
            return self._create_project
        if method == "POST" and route == "/project":
            return self._create_project

        # 实验
        if method == "POST" and len(parts) == 4 and parts[0] == "project" and parts[3] == "experiment":
            return self._create_experiment
        if method == "GET" and len(parts) == 5 and parts[0] == "project" and parts[3] == "runs":
            return self._get_run
        if method == "PUT" and len(parts) == 6 and parts[0] == "project" and parts[3] == "runs" and parts[5] == "state":
            return self._stop_experiment
        if method == "PUT" and len(parts) == 6 and parts[0] == "project" and parts[5] == "profile":
            return self._profile
        if method == "PUT" and len(parts) == 6 and parts[0] == "projects" and parts[5] == "profile":
            return self._profile

        # 指标
        if method == "GET" and len(parts) == 5 and parts[:2] == ["house", "metrics"] and parts[2] == "summaries":
            return self._summary
        if method == "POST" and route == "/house/metrics":
            return self._metrics
        if method == "POST" and len(parts) == 4 and parts[0] == "house" and parts[3] == "heartbeat":
            return self._heartbeat
        if method == "POST" and len(parts) == 4 and parts[0] == "projects" and parts[3] == "series":
            return self._columns

        # 文件
        if method == "POST" and route == "/resources/presigned/put":
            return self._presigned
        if method == "POST" and len(parts) == 4 and parts[0] == "experiment" and parts[2] == "files":
            return getattr(self, "_files_" + parts[3].replace("-", "_"), None)
        return None

    # ------------------------------------------------------------ 处理函数

    def _verify(self, _body: bytes, _route: str) -> None:
        self._json(200, VERIFY_PROFILE)

    def _health(self, _body: bytes, _route: str) -> None:
        self._json(200, {"status": "ok"})

    def _login(self, _body: bytes, _route: str) -> None:
        self._json(
            200,
            {
                "sid": "mock-sid",
                "expiredAt": "2999-01-01T00:00:00.000000Z",
                "userInfo": {"username": "bench"},
            },
        )

    def _get_project(self, _body: bytes, route: str) -> None:
        parts = [p for p in route.split("/") if p]
        owner, name = parts[1], parts[2]
        if name == "brand-new" and name not in self.state.projects:
            self._json(404, {"code": "not_found", "message": "project not found"})
            return
        with self.state.lock:
            experiments_seen = self.state.experiment_seq
        self._json(
            200,
            {
                "cuid": f"proj-{name}",
                "name": name,
                "path": f"/{owner}/{name}",
                "visibility": "PRIVATE",
                "version": 1,
                "group": {"username": owner},
                # 实验名与配色由后端实验数量派生，这里给一个固定值。
                "_count": {"experiments": experiments_seen, "contributors": 1, "collaborators": 0, "clones": 0},
            },
        )

    def _create_project(self, body: bytes, route: str) -> None:
        payload = _decode(body)
        name = payload.get("name") or route.rsplit("/", 1)[-1]
        with self.state.lock:
            self.state.projects.add(name)
        self._json(201, {"cuid": f"proj-{name}"})

    def _create_experiment(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        run_id = str(payload.get("cuid") or payload.get("id") or "")
        if not run_id:
            run_id = f"anon-{int(time.monotonic() * 1000)}"
        experiment_id, created = self.state.ensure_experiment(run_id)
        if created:
            # 新建实验时后端不回传 createdAt，SDK 会用上报时刻补齐。
            self._json(201, {"cuid": experiment_id, "slug": run_id, "name": payload.get("name", "")})
            return
        self._json(
            200,
            {
                "cuid": experiment_id,
                "slug": run_id,
                "name": payload.get("name", ""),
                "createdAt": "2026-01-02T00:00:00.000Z",
            },
        )

    def _get_run(self, _body: bytes, route: str) -> None:
        parts = [p for p in route.split("/") if p]
        run_id = parts[4]
        self._json(200, {"cuid": run_id, "createdAt": "2026-01-03T00:00:00.000Z"})

    def _stop_experiment(self, _body: bytes, _route: str) -> None:
        with self.state.lock:
            self.state.stop_experiments += 1
        self._json(200, {"ok": True})

    def _profile(self, _body: bytes, _route: str) -> None:
        with self.state.lock:
            self.state.profile_uploads += 1
        self._json(200, {"ok": True})

    def _summary(self, _body: bytes, _route: str) -> None:
        # 空摘要：非新建 start 会拉取，返回空表示云端还没有任何指标。
        self._json(200, {"scalar": [], "media": [], "log": [], "histogram": []})

    def _metrics(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        metric_type = str(payload.get("type") or "unknown")
        metrics = payload.get("metrics") or []
        if not isinstance(metrics, list):
            metrics = []
        self.state.ingest_metrics(metric_type, metrics, len(body), 200)
        self._json(200, {"ok": True})

    def _heartbeat(self, _body: bytes, _route: str) -> None:
        with self.state.lock:
            self.state.heartbeats += 1
        self._json(200, {"ok": True})

    def _columns(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        series = payload.get("series") or []
        keys = [str(item.get("key")) for item in series if isinstance(item, dict)]
        self.state.ingest_columns(keys)
        self._json(200, {"ok": True})

    def _presigned(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        paths = payload.get("paths") or []
        port = self._port()
        with self.state.lock:
            self.state.presigned_requests += 1
            self.state.presigned_paths += len(paths)
        self._json(200, {"urls": [f"http://127.0.0.1:{port}/sink/{path}" for path in paths]})

    def _files_prepare(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        files = payload.get("files") or []
        port = self._port()
        with self.state.lock:
            self.state.files_prepared += len(files)
        urls = [f"http://127.0.0.1:{port}/sink/{item.get('path') or f'save-{i}'}" for i, item in enumerate(files)]
        self._json(200, {"urls": urls})

    def _files_complete(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        files = payload.get("files") or []
        with self.state.lock:
            self.state.files_completed += len(files)
        self._json(201, {"ok": True})

    def _files_prepare_multipart(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        files = payload.get("files") or []
        port = self._port()
        with self.state.lock:
            self.state.files_prepared += len(files)
        out = []
        for index, item in enumerate(files):
            count = int(item.get("count") or 1)
            base = item.get("path") or f"multipart-{index}"
            out.append(
                {
                    "uploadId": f"upload-{index}",
                    "parts": [
                        {"partNumber": part, "url": f"http://127.0.0.1:{port}/sink/{base}.part{part}"}
                        for part in range(1, count + 1)
                    ],
                }
            )
        self._json(200, {"files": out})

    def _files_complete_multipart(self, body: bytes, _route: str) -> None:
        payload = _decode(body)
        files = payload.get("files") or []
        with self.state.lock:
            self.state.files_completed += len(files)
        out = [{"path": item.get("path", ""), "state": "UPLOADED"} for item in files]
        self._json(201, {"files": out})


def _decode(body: bytes) -> Dict[str, Any]:
    if not body:
        return {}
    try:
        payload = json.loads(body)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


class MockServer(ThreadingHTTPServer):
    """默认 backlog 只有 5，多 rank 并发上传会被重置。"""

    daemon_threads = True
    allow_reuse_address = True
    request_queue_size = 256


def serve(port: int, port_file: Optional[str], latency_ms: int) -> int:
    handler = type("BoundMockHandler", (MockHandler,), {"state": CloudState(), "latency_ms": latency_ms})
    server = MockServer(("127.0.0.1", port), handler)
    actual_port = server.server_address[1]
    if port_file:
        tmp = f"{port_file}.tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            handle.write(str(actual_port))
        os.replace(tmp, port_file)
    print(f"MOCK_PORT={actual_port}", flush=True)
    try:
        server.serve_forever(poll_interval=0.02)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(description="SwanLab bench mock cloud (counting only)")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--port-file", default=None)
    parser.add_argument("--latency-ms", type=int, default=0)
    args = parser.parse_args(argv)
    return serve(args.port, args.port_file, args.latency_ms)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
