"""视频本地 HTTP 回归：流式 PUT、逐文件结果、JSON 重试和 Transport 隔离。"""

import hashlib
import json
import threading
import time
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from typing import Literal
from urllib.parse import urlsplit

import pytest
import requests
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab import Video
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import (
    HistogramRecord,
    MediaItem,
    MediaRecord,
    MediaValue,
    ScalarRecord,
)
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.transport.sender.media import collect_media_tasks
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.transport.worker import Transport
from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy


class Backend:
    def __init__(self):
        self.metrics = []
        self.signatures = []
        self.puts = []
        self.columns = []
        self.gets = []
        self.failures = {}
        self.metric_failures = 0
        self.metric_status = 200
        self.started = threading.Event()
        self.release = threading.Event()
        self.release.set()
        self.active = 0
        self.max_active = 0
        backend = self

        class Handler(BaseHTTPRequestHandler):
            def reply(self, value, status=200):
                content = json.dumps(value).encode()
                self.send_response(status)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(content)))
                self.end_headers()
                self.wfile.write(content)

            def do_POST(self):
                data = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
                if self.path.endswith("/login/api_key"):
                    expired = (datetime.now(timezone.utc) + timedelta(days=30)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                    self.reply({"sid": "test-sid", "expiredAt": expired, "userInfo": {"username": "tester"}})
                elif self.path.endswith("/resources/presigned/put"):
                    backend.signatures.append(data)
                    self.reply({"urls": [f"{backend.url}/objects/{path}?signature=test" for path in data["paths"]]})
                elif self.path.endswith("/series"):
                    backend.columns.append(data)
                    self.reply({}, 201)
                elif self.path.endswith("/house/metrics"):
                    backend.metrics.append(data)
                    status = 500 if data["type"] == "media" and backend.metric_failures else backend.metric_status
                    if status == 500:
                        backend.metric_failures -= 1
                    self.reply("SUCCESS" if status == 200 else {"code": "unavailable"}, status)
                else:
                    self.reply({}, 404)

            def do_PUT(self):
                path = urlsplit(self.path).path.removeprefix("/objects/")
                is_video = "/video/" in path
                if is_video:
                    backend.active += 1
                    backend.max_active = max(backend.max_active, backend.active)
                    backend.started.set()
                    backend.release.wait(timeout=5)
                try:
                    digest = hashlib.sha256()
                    remaining = int(self.headers["Content-Length"])
                    size = remaining
                    while remaining:
                        chunk = self.rfile.read(min(65536, remaining))
                        if not chunk:
                            break
                        digest.update(chunk)
                        remaining -= len(chunk)
                    backend.puts.append((path, self.headers["Content-Type"], size, digest.hexdigest()))
                    statuses = backend.failures.get(path, [])
                    status = statuses.pop(0) if statuses else 200
                    self.reply({"code": "ExpiredToken"} if status == 401 else {}, status)
                finally:
                    if is_video:
                        backend.active -= 1

            def do_GET(self):
                backend.gets.append(self.path)
                self.reply({}, 404)

            def log_message(self, format, *args):
                pass

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.server.server_port}"
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def close(self):
        self.release.set()
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)


@pytest.fixture
def backend(monkeypatch):
    server = Backend()
    if client.exists():
        client.reset()
    client.new(api_key="test-key", base_url=server.url)
    monkeypatch.setattr(RetryPolicy, "schedule", lambda self, rounds: 0.005)
    try:
        yield server
    finally:
        client.reset()
        server.close()


def context(path, mode: Literal["core", "sync"] = "core", *, record_batch=1000, skip_store=False):
    ctx = CoreContext(
        config=CoreConfig(
            run_id="video-run",
            run_dir=path,
            section_rule=0,
            record_batch=record_batch,
            record_interval=0.01,
            save_size=1024**3,
            save_split=1024**2,
            save_part=1024**2,
            save_batch=100,
            skip_store=skip_store,
        ),
        mode=mode,
    )
    ctx.set_online_params("alice", "project", "project-id", 1, "experiment-id")
    return ctx


def video_record(ctx, count=1, key="clips", step=0, number=1):
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    items = []
    for index in range(count):
        item = Video(b"GIF89a" + bytes([index + 1]), caption="" if index == 0 else str(index)).transform(
            step=step,
            path=None if ctx.config.skip_store else ctx.media_dir / "video",
        )
        assert item is not None
        items.append(item)
    return Record(
        num=number,
        media=MediaRecord(
            key=key, step=step, type=ColumnType.COLUMN_TYPE_VIDEO, timestamp=timestamp, value=MediaValue(items=items)
        ),
    )


def drain(sender, records):
    pending = list(records)
    deadline = time.monotonic() + 5
    while pending and time.monotonic() < deadline:
        result = sender.upload_media(pending)
        pending = result.remaining_records
        time.sleep(0.005)
    assert not pending


def test_partial_put_failures_and_stable_json_retry(tmp_path, backend):
    ctx = context(tmp_path)
    record = video_record(ctx, count=3)
    paths = [f"media/video/{item.filename}" for item in record.media.value.items]
    backend.failures[paths[1]] = [403]
    backend.failures[paths[2]] = [500, 200]
    backend.metric_failures = 1
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert [put[0] for put in backend.puts].count(paths[0]) == 1
        assert [put[0] for put in backend.puts].count(paths[1]) == 1
        assert [put[0] for put in backend.puts].count(paths[2]) == 2
        assert len(backend.metrics) == 2
        assert backend.metrics[0] == backend.metrics[1]
        metric = backend.metrics[0]["metrics"][0]
        assert metric["data"] == [paths[0], paths[2]]
        assert metric["more"] == [{}, {"caption": "2"}]
        assert all(len(request["paths"]) == 1 for request in backend.signatures)
        assert backend.max_active == 1 and not backend.gets
        assert all(put[1] == "image/gif" for put in backend.puts)
    finally:
        sender.close()


def test_expired_signature_reacquired(tmp_path, backend):
    ctx = context(tmp_path)
    record = video_record(ctx)
    path = f"media/video/{record.media.value.items[0].filename}"
    backend.failures[path] = [401, 200]
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert len(backend.signatures) == len(backend.puts) == 2
        assert len(backend.metrics) == 1
    finally:
        sender.close()


def test_video_above_ordinary_media_limit_streams_unchanged(tmp_path, backend):
    ctx = context(tmp_path)
    content = b"GIF89a" + b"x" * (11 * 1024**2)
    item = Video(content).transform(step=0, path=ctx.media_dir / "video")
    assert item is not None
    ts = Timestamp()
    ts.GetCurrentTime()
    record = Record(
        num=1,
        media=MediaRecord(key="large", type=ColumnType.COLUMN_TYPE_VIDEO, timestamp=ts, value=MediaValue(items=[item])),
    )
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert backend.puts[0][2:] == (len(content), hashlib.sha256(content).hexdigest())
        assert backend.metrics[0]["metrics"][0]["more"] == [{}]
    finally:
        sender.close()


@pytest.mark.parametrize("mode", ["core", "sync"])
def test_blocked_video_put_does_not_block_other_metrics(tmp_path, backend, mode):
    ctx = context(tmp_path, mode)
    backend.release.clear()
    tracker = UploadTracker()
    transport = Transport(ctx, tracker=tracker)
    record = video_record(ctx)
    try:
        transport.put(
            [
                Record(num=2, column=ColumnRecord(column_key="clips", column_type=ColumnType.COLUMN_TYPE_VIDEO)),
                record,
                video_record(ctx, key="other", step=1, number=3),
            ]
        )
        assert backend.started.wait(timeout=3)
        ts = Timestamp()
        ts.GetCurrentTime()
        transport.put(
            [
                Record(num=4, scalar=ScalarRecord(key="loss", step=0, timestamp=ts, value={"number": 0.5})),
                Record(
                    num=5,
                    histogram=HistogramRecord(
                        key="hist", step=0, timestamp=ts, value={"values": [1], "bins": [0.0, 1.0]}
                    ),
                ),
                Record(
                    num=6,
                    media=MediaRecord(
                        key="image",
                        type=ColumnType.COLUMN_TYPE_IMAGE,
                        timestamp=ts,
                        value=MediaValue(items=[MediaItem(filename="test.gif", payload=b"GIF89a", size=6)]),
                    ),
                ),
            ]
        )
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline and len(backend.metrics) < 3:
            time.sleep(0.01)
        assert {metric["type"] for metric in backend.metrics} == {"scalar", "histogram", "media"}
        assert transport.pending_records()
        backend.release.set()
        assert transport.finish(timeout=5), transport.finish_summary()
        assert backend.max_active == 1
        assert not transport.pending_records()
        stats = tracker.snapshot()
        assert stats.total_records == stats.uploaded_records == 6
        assert backend.columns[0]["series"][0]["type"] == "VIDEO"
    finally:
        backend.release.set()
        transport.finish(timeout=5)


def test_skip_store_http_retry_sends_identical_full_body(tmp_path, backend):
    ctx = context(tmp_path, skip_store=True)
    record = video_record(ctx)
    item = record.media.value.items[0]
    backend.failures[f"media/video/{item.filename}"] = [500, 200]
    tracker = UploadTracker()
    sender = HttpRecordSender(ctx)
    sender.set_tracker(tracker)
    try:
        drain(sender, [record])
        expected = (len(item.payload), hashlib.sha256(item.payload).hexdigest())
        assert [put[2:] for put in backend.puts] == [expected, expected]
        assert not tracker._files
        assert not ctx.media_dir.exists()
    finally:
        sender.close()


def test_skip_store_partial_read_failure_rewinds_before_retry(tmp_path, monkeypatch):
    from swanlab.sdk.internal.core_python.transport.sender import video as sender_video

    ctx = context(tmp_path, skip_store=True)
    record = video_record(ctx)
    uploader = sender_video.VideoUploader(ctx, UploadTracker())
    reads = []

    def put(session, experiment, payloads, **kwargs):
        payload = payloads[0]
        assert isinstance(payload.buffer, BytesIO)
        if not reads:
            reads.append(payload.buffer.read(3))
            return {0: requests.ConnectionError("interrupted")}
        reads.append(payload.buffer.read())
        return {}

    monkeypatch.setattr(sender_video, "put_media_resources", put)
    monkeypatch.setattr(sender_video, "upload_media", lambda *args, **kwargs: None)
    uploader._retry.schedule = lambda failed_rounds: 0.0
    try:
        assert uploader.upload([record]).remaining_records
        assert uploader.upload([record]).confirmed_records
        assert reads == [record.media.value.items[0].payload[:3], record.media.value.items[0].payload]
    finally:
        uploader.close()


def test_backoff_in_first_record_does_not_block_later_dispatch_chunks(tmp_path, backend):
    ctx = context(tmp_path, record_batch=1)
    first = video_record(ctx, number=1)
    second = video_record(ctx, key="second", number=2, step=2)
    first_path = f"media/video/{first.media.value.items[0].filename}"
    second_path = f"media/video/{second.media.value.items[0].filename}"
    backend.failures[first_path] = [500]
    sender = HttpRecordSender(ctx)
    dispatch = Dispatch(batch_size=1, sender=sender)
    try:
        result = dispatch([first, second])
        assert len(backend.puts) == 1
        assert sender._video_uploader is not None
        state = sender._video_uploader._tasks[id(first)]
        state.retry_at[0] = time.monotonic() + 30
        result = dispatch(result.remaining_records)
        assert [put[0] for put in backend.puts] == [first_path, second_path]
        assert result.remaining_records == [first]
        assert result.confirmed_records == [second]
    finally:
        sender.close()


def test_new_video_record_runs_while_earlier_record_is_in_backoff(tmp_path, backend, monkeypatch):
    monkeypatch.setattr(RetryPolicy, "schedule", lambda self, rounds: 5)
    ctx = context(tmp_path, record_batch=1)
    first = video_record(ctx)
    backend.failures[f"media/video/{first.media.value.items[0].filename}"] = [500, 200]
    transport = Transport(ctx)
    try:
        transport.put([first])
        assert backend.started.wait(timeout=3)
        deadline = time.monotonic() + 3
        while not backend.puts and time.monotonic() < deadline:
            time.sleep(0.005)
        transport.put([video_record(ctx, key="second", step=2, number=2)])
        deadline = time.monotonic() + 2
        while len(backend.puts) < 2 and time.monotonic() < deadline:
            time.sleep(0.005)
        assert len(backend.puts) >= 2
        assert backend.puts[1][0] != backend.puts[0][0]
        assert transport._video is not None and transport._video._sender is not None
        uploader = transport._video._sender._video_uploader
        assert uploader is not None
        uploader._tasks[id(first)].retry_at[0] = 0
        assert transport.finish(timeout=5)
    finally:
        if transport._video is not None and transport._video._sender is not None:
            uploader = transport._video._sender._video_uploader
            if uploader is not None:
                for task in uploader._tasks.values():
                    task.retry_at.clear()
        transport.finish(timeout=5)


def test_permanent_metric_rejection_closes_tracker_and_finishes(tmp_path, backend):
    ctx = context(tmp_path)
    backend.metric_status = 403
    tracker = UploadTracker()
    transport = Transport(ctx, tracker=tracker)
    record = video_record(ctx)
    try:
        transport.put([record])
        assert transport.finish(timeout=5)
        assert not transport.pending_records()
        assert not tracker._files
        assert (ctx.media_dir / "video" / record.media.value.items[0].filename).exists()
    finally:
        transport.finish(timeout=5)


@pytest.mark.parametrize(
    "media_type", [ColumnType.COLUMN_TYPE_VIDEO, ColumnType.COLUMN_TYPE_IMAGE, ColumnType.COLUMN_TYPE_AUDIO]
)
def test_caption_alignment_is_limited_to_video(tmp_path, media_type):
    ctx = context(tmp_path, skip_store=True)
    ts = Timestamp()
    ts.GetCurrentTime()
    record = Record(
        media=MediaRecord(
            type=media_type,
            timestamp=ts,
            value=MediaValue(
                items=[
                    MediaItem(filename="first.gif", payload=b"GIF89a", size=6),
                    MediaItem(filename="second.gif", payload=b"GIF89a", size=6, caption="second"),
                ]
            ),
        )
    )
    task = collect_media_tasks(ctx, [record])[0]
    assert task.metric is not None
    expected = [{}, {"caption": "second"}] if media_type == ColumnType.COLUMN_TYPE_VIDEO else [{"caption": "second"}]
    assert task.metric["more"] == expected
