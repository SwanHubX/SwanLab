"""视频本地 HTTP 回归：复用媒体上传、record 重试与后台线程清理。"""

import hashlib
import json
import threading
import time
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from typing import Literal
from urllib.parse import urlsplit
from uuid import uuid4

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
from swanlab.sdk.internal.core_python.transport import sender as sender_module
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.transport.sender import media as sender_media
from swanlab.sdk.internal.core_python.transport.sender.media import collect_media_tasks
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.transport.worker import Transport
from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy
from swanlab.sdk.internal.pkg import adapter


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
        self.lock = threading.Lock()
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
                    with backend.lock:
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
                        with backend.lock:
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
        assert not self.thread.is_alive()

    def metric_requests(self, ctx):
        return [
            request
            for request in self.metrics
            if request["projectId"] == ctx.project_id and request["experimentId"] == ctx.experiment_id
        ]

    def media_items(self, ctx, key):
        return [
            metric
            for request in self.metric_requests(ctx)
            if request["type"] == "media"
            for metric in request["metrics"]
            if metric["key"] == key
        ]


@pytest.fixture
def backend(monkeypatch, transport_cleanup):
    server = Backend()
    if client.exists():
        client.reset()
    client.new(api_key="test-key", base_url=server.url)
    monkeypatch.setattr(RetryPolicy, "schedule", lambda self, rounds: 0.005)
    try:
        yield server
    finally:
        server.release.set()
        try:
            transport_cleanup()
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
    ctx.set_online_params("alice", "project", "video-project", 1, uuid4().hex[:21])
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


def test_media_record_retry_reuploads_files_and_preserves_metadata(tmp_path, backend):
    ctx = context(tmp_path)
    record = video_record(ctx, count=3)
    paths = [f"media/video/{item.filename}" for item in record.media.value.items]
    backend.failures[paths[1]] = [403]
    backend.failures[paths[2]] = [500, 200]
    backend.metric_failures = 1
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert all([put[0] for put in backend.puts].count(path) == 3 for path in paths)
        requests = backend.metric_requests(ctx)
        assert len(requests) == 2
        assert requests[0] == requests[1]
        metric = backend.media_items(ctx, "clips")[0]
        assert metric["data"] == paths
        assert metric["more"] == [{}, {"caption": "1"}, {"caption": "2"}]
        assert all(request["paths"] == paths for request in backend.signatures)
        assert backend.max_active <= 3 and not backend.gets
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
        assert len(backend.metric_requests(ctx)) == 1
    finally:
        sender.close()


@pytest.mark.parametrize("skip_store", [False, True])
@pytest.mark.parametrize(
    "fmt,content,mime",
    [
        ("gif", b"GIF89a-payload", "image/gif"),
        ("mp4", b"\x00\x00\x00\x10ftypmp42\x00\x00\x00\x00", "video/mp4"),
        ("webm", b"\x1a\x45\xdf\xa3\x87\x42\x82\x84webm", "video/webm"),
        ("ogg", b"OggS-payload", "video/ogg"),
    ],
)
def test_video_formats_use_media_upload(tmp_path, backend, skip_store, fmt, content, mime):
    ctx = context(tmp_path, skip_store=skip_store)
    item = Video(content, format=fmt).transform(step=0, path=None if skip_store else ctx.media_dir / "video")
    assert item is not None
    ts = Timestamp()
    ts.GetCurrentTime()
    record = Record(
        media=MediaRecord(
            key="formats", type=ColumnType.COLUMN_TYPE_VIDEO, timestamp=ts, value=MediaValue(items=[item])
        )
    )
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert backend.puts == [
            (f"media/video/{item.filename}", mime, len(content), hashlib.sha256(content).hexdigest())
        ]
        assert backend.media_items(ctx, "formats")[0]["more"] == [{}]
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
        metrics = backend.media_items(ctx, "large")
        assert len(metrics) == 1
        assert metrics[0]["more"] == [{}]
    finally:
        sender.close()


def test_video_metadata_is_selected_by_experiment_and_key(tmp_path, backend):
    ctx = context(tmp_path)
    for experiment, metric_type, key in [
        ("other-experiment", "media", "clips"),
        (ctx.experiment_id, "scalar", "clips"),
        (ctx.experiment_id, "media", "other-key"),
    ]:
        client.post(
            "/house/metrics",
            {
                "projectId": ctx.project_id,
                "experimentId": experiment,
                "type": metric_type,
                "metrics": [{"key": key, "data": 1}],
            },
            retries=0,
        )
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [video_record(ctx)])
        metrics = backend.media_items(ctx, "clips")
        assert len(metrics) == 1
        assert metrics[0]["more"] == [{}]
    finally:
        sender.close()


@pytest.mark.parametrize("mode", ["core", "sync"])
def test_video_and_other_metrics_share_upload_thread(tmp_path, backend, mode):
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
        deadline = time.monotonic() + 3
        while backend.active < 2 and time.monotonic() < deadline:
            time.sleep(0.005)
        assert backend.active == 2
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
        assert not backend.metric_requests(ctx)
        assert transport._thread is not None and transport._thread.is_alive()
        assert sum(thread.name == Transport.THREAD_NAME for thread in threading.enumerate()) == 1
        assert {record.num for record in transport.pending_records()} >= {1, 3, 4, 5, 6}
        assert not transport.finish(timeout=0.01)
        backend.release.set()
        assert transport.finish(timeout=5), transport.finish_summary()
        assert {metric["type"] for metric in backend.metric_requests(ctx)} == {"scalar", "histogram", "media"}
        assert backend.max_active == 2
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
    ctx = context(tmp_path, skip_store=True)
    record = video_record(ctx)
    sender = HttpRecordSender(ctx)
    reads = []

    def put(session, experiment, *, buffers, **kwargs):
        buffer = buffers[0]
        assert isinstance(buffer, BytesIO)
        if not reads:
            reads.append(buffer.read(3))
            return {0: requests.ConnectionError("interrupted")}
        reads.append(buffer.read())
        return {}

    monkeypatch.setattr(sender_media, "upload_resource", put)
    monkeypatch.setattr(sender_module, "upload_media", lambda *args, **kwargs: None)
    try:
        assert sender.upload_media([record]).remaining_records
        assert sender.upload_media([record]).confirmed_records
        assert reads == [record.media.value.items[0].payload[:3], record.media.value.items[0].payload]
    finally:
        sender.close()


def test_video_record_retry_respects_media_dispatch_batches(tmp_path, backend):
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
        assert result.remaining_records == [first, second]
        result = dispatch(result.remaining_records)
        assert [put[0] for put in backend.puts] == [first_path, first_path, second_path]
        assert not result.remaining_records
        assert result.confirmed_records == [first, second]
    finally:
        sender.close()


def test_new_video_record_waits_for_pending_upload(tmp_path, backend, monkeypatch):
    monkeypatch.setattr(RetryPolicy, "schedule", lambda self, rounds: 0.2)
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
        assert transport.finish(timeout=5)
        assert len(backend.puts) == 3
        assert backend.puts[0][0] == backend.puts[1][0]
        assert backend.puts[2][0] != backend.puts[0][0]
    finally:
        assert transport.finish(timeout=5)


@pytest.mark.parametrize("skip_store", [False, True])
def test_permanent_metric_rejection_closes_tracker_and_finishes(tmp_path, backend, skip_store, monkeypatch):
    from swanlab.sdk.internal.pkg import console

    warnings = []
    monkeypatch.setattr(console, "warning", lambda message: warnings.append(message))
    ctx = context(tmp_path, skip_store=skip_store)
    backend.metric_status = 403
    tracker = UploadTracker()
    transport = Transport(ctx, tracker=tracker)
    record = video_record(ctx)
    try:
        transport.put([record])
        assert transport.finish(timeout=5)
        assert not transport.pending_records()
        assert not tracker._files
        if skip_store:
            assert not ctx.media_dir.exists()
            assert not any("local video files retained" in message for message in warnings)
        else:
            assert (ctx.media_dir / "video" / record.media.value.items[0].filename).exists()
    finally:
        transport.finish(timeout=5)


@pytest.mark.parametrize(
    "media_type",
    [
        ColumnType.COLUMN_TYPE_VIDEO,
        ColumnType.COLUMN_TYPE_IMAGE,
        ColumnType.COLUMN_TYPE_AUDIO,
        ColumnType.COLUMN_TYPE_TEXT,
        ColumnType.COLUMN_TYPE_HTML,
        ColumnType.COLUMN_TYPE_ECHARTS,
        ColumnType.COLUMN_TYPE_OBJECT3D,
        ColumnType.COLUMN_TYPE_MOLECULE,
    ],
)
@pytest.mark.parametrize("skip_store", [False, True])
def test_caption_alignment_for_all_media(tmp_path, media_type, skip_store):
    ctx = context(tmp_path, skip_store=skip_store)
    ts = Timestamp()
    ts.GetCurrentTime()
    items = [
        MediaItem(filename="first.gif", size=6),
        MediaItem(filename="missing.gif", size=6, caption="missing"),
        MediaItem(filename="second.gif", size=6, caption="second"),
    ]
    if skip_store:
        items[0].payload = items[2].payload = b"GIF89a"
    record = Record(
        media=MediaRecord(
            type=media_type,
            timestamp=ts,
            value=MediaValue(items=items),
        )
    )
    if not skip_store:
        path = ctx.media_dir / adapter.medium[media_type]
        path.mkdir(parents=True)
        (path / "first.gif").write_bytes(b"GIF89a")
        (path / "second.gif").write_bytes(b"GIF89a")
    task = collect_media_tasks(ctx, [record])[0]
    assert task.metric is not None
    assert task.metric["more"] == [{}, {"caption": "second"}]
    assert task.metric["data"] == [f"media/{adapter.medium[media_type]}/{name}.gif" for name in ("first", "second")]


def test_image_caption_alignment_in_actual_http_request(tmp_path, backend):
    ctx = context(tmp_path, skip_store=True)
    ts = Timestamp()
    ts.GetCurrentTime()
    record = Record(
        media=MediaRecord(
            key="images",
            type=ColumnType.COLUMN_TYPE_IMAGE,
            timestamp=ts,
            value=MediaValue(
                items=[
                    MediaItem(filename="first.png", payload=b"first", size=5),
                    MediaItem(filename="second.png", payload=b"second", size=6, caption="b"),
                ]
            ),
        )
    )
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        metrics = backend.media_items(ctx, "images")
        assert len(metrics) == 1
        assert metrics[0]["data"] == ["media/image/first.png", "media/image/second.png"]
        assert metrics[0]["more"] == [{}, {"caption": "b"}]
    finally:
        sender.close()


def test_video_record_retries_do_not_rehash_in_sender(tmp_path, backend, monkeypatch):
    ctx = context(tmp_path)
    record = video_record(ctx)
    backend.failures[f"media/video/{record.media.value.items[0].filename}"] = [500, 200]
    backend.metric_failures = 1
    original = hashlib.sha256
    sender_hashes = []

    def sha256(*args, **kwargs):
        if threading.current_thread() is threading.main_thread():
            sender_hashes.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(hashlib, "sha256", sha256)
    sender = HttpRecordSender(ctx)
    try:
        drain(sender, [record])
        assert len(backend.puts) == 3
        assert len(backend.media_items(ctx, "clips")) == 2
        assert not sender_hashes
    finally:
        sender.close()
