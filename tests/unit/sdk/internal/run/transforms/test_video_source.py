"""视频快照、大小边界、共享 buffer 和错误隔离。"""

import gc
import hashlib
import inspect
from io import BytesIO
from types import SimpleNamespace
from typing import cast
from unittest.mock import Mock

import numpy as np
import pytest
from google.protobuf.timestamp_pb2 import Timestamp

import swanlab
import swanlab.sdk.internal.run.transforms.video as video_module
from swanlab import Video
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord, MediaValue
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run import Run
from swanlab.sdk.internal.run.components.consumer.builder import RecordBuilder
from swanlab.sdk.internal.run.transforms.video import source

GIF = b"GIF89a" + b"payload"


def test_bytesio_cursor_and_invalid_offset():
    stream = BytesIO(b"prefix" + GIF)
    stream.seek(6)
    video = Video(stream)
    assert video.is_valid and stream.tell() == 6 + len(GIF)
    assert video.buffer.getvalue() == GIF
    invalid = BytesIO(GIF)
    invalid.seek(2)
    assert not Video(invalid).is_valid
    assert invalid.tell() == len(GIF)


def test_binding_errors_remain_python_errors():
    with pytest.raises(TypeError):
        Video()  # type: ignore
    with pytest.raises(TypeError):
        Video(GIF, unknown=True)  # type: ignore


def test_snapshot_streaming_lazy_buffer_and_source_preserved(tmp_path, monkeypatch):
    original = tmp_path / "user.gif"
    content = GIF + b"a" * (source.CHUNK_SIZE * 2 + 17)
    original.write_bytes(content)
    scratch = tmp_path / ".scratch" / "owned.gif"
    scratch.parent.mkdir()
    scratch.touch()
    monkeypatch.setattr(video_module, "scratch_path", lambda _: scratch)
    video = Video(original)
    assert video.is_valid and video._source is not None
    assert video._source.memory is None
    original.write_bytes(GIF + b"changed")
    monkeypatch.setattr(
        source.VideoSource, "buffer", property(lambda _: pytest.fail("default path must not access buffer"))
    )
    item = video.transform(step=2, path=tmp_path / "media")
    assert item is not None
    assert item.sha256 == hashlib.sha256(content).hexdigest()
    assert (tmp_path / "media" / item.filename).read_bytes() == content
    assert not scratch.exists()
    assert original.read_bytes() == GIF + b"changed"


def test_nested_source_shared_buffer_changes_rehashed(tmp_path):
    inner = Video(GIF, "caption")
    outer = Video(inner)
    assert outer.buffer is inner.buffer
    outer.buffer.seek(0, 2)
    outer.buffer.write(b"changed")
    item = inner.transform(step=3, path=tmp_path)
    assert item is not None
    assert item.sha256 == hashlib.sha256(GIF + b"changed").hexdigest()
    assert outer.caption == "caption"
    assert inner.sha256 == outer.sha256


def test_owned_scratch_released_with_last_nested_source(tmp_path, monkeypatch):
    scratch = tmp_path / "owned.gif"
    scratch.touch()
    monkeypatch.setattr(video_module, "scratch_path", lambda _: scratch)
    inner = Video(GIF)
    outer = Video(inner)
    del inner
    gc.collect()
    assert scratch.exists()
    del outer
    gc.collect()
    assert not scratch.exists()


@pytest.mark.parametrize("delta,valid", [(-1, True), (0, True), (1, False)])
def test_exact_size_limit(delta, valid, monkeypatch):
    monkeypatch.setattr(source, "MAX_SIZE", 100)
    content = GIF + b"x" * (100 + delta - len(GIF))
    video = Video(content)
    assert video.is_valid is valid


def test_buffer_growth_rechecked_before_logging(monkeypatch):
    monkeypatch.setattr(source, "MAX_SIZE", 100)
    video = Video(GIF)
    video.buffer.seek(0, 2)
    video.buffer.write(b"x" * 100)
    assert video.transform(step=0, path=None) is None
    assert video.error is not None and video.error.code == "size_limit"


def test_media_limits_remain_separate(tmp_path):
    builder = RecordBuilder(cast(RunContext, SimpleNamespace(media_dir=tmp_path)))
    video = MediaRecord(
        key="video", type=ColumnType.COLUMN_TYPE_VIDEO, value=MediaValue(items=[MediaItem(size=11 * 1024**2)])
    )
    image = MediaRecord(key="image", type=ColumnType.COLUMN_TYPE_IMAGE, value=video.value)
    assert builder._ensure_media_size(video) is not None
    assert builder._ensure_media_size(image) is None
    video.value.items[0].size = source.MAX_SIZE
    assert builder._ensure_media_size(video) is not None
    video.value.items[0].size += 1
    assert builder._ensure_media_size(video) is None


def test_invalid_items_filtered_before_108_limit(tmp_path):
    ctx = SimpleNamespace(
        media_dir=tmp_path, config=SimpleNamespace(settings=SimpleNamespace(core=SimpleNamespace(skip_store=True)))
    )
    builder = RecordBuilder(cast(RunContext, ctx))
    values = [Video(b"invalid")] + [Video(GIF, caption=str(i)) for i in range(109)]
    record, _ = builder.build_scalar_or_media(values, "clips", Timestamp(), 0)
    assert isinstance(record, MediaRecord)
    assert len(record.value.items) == 108
    assert record.value.items[0].caption == "0"
    assert record.value.items[-1].caption == "107"
    empty, _ = builder.build_scalar_or_media([Video(b"invalid")], "empty", Timestamp(), 0)
    assert empty is None


def test_log_video_keeps_5d_array_as_one_video():
    run = SimpleNamespace(log=Mock())
    frames = np.zeros((3, 2, 3, 4, 4), dtype=np.uint8)
    inspect.unwrap(Run.log_video)(run, key="batch", data=frames, fps=8, format="gif", step=7)
    payload = run.log.call_args.args[0]
    assert len(payload["batch"]) == 1
    assert payload["batch"][0].is_valid
    assert run.log.call_args.kwargs["step"] == 7


def test_log_video_caption_length_failure_is_local_to_call():
    run = SimpleNamespace(log=Mock())
    inspect.unwrap(Run.log_video)(run, key="clips", data=[GIF], caption=["one", "two"])
    run.log.assert_not_called()
    inspect.unwrap(Run.log_video)(run, key="good", data=GIF)
    run.log.assert_called_once()


@pytest.mark.parametrize("input_kind", ["path", "bytes", "array"])
def test_disabled_video_has_no_sdk_file_side_effects(tmp_path, input_kind):
    path = tmp_path / "user.gif"
    path.write_bytes(GIF)
    data = path if input_kind == "path" else GIF if input_kind == "bytes" else np.zeros((2, 3, 8, 8), dtype=np.uint8)
    run = swanlab.init(mode="disabled", logdir=str(tmp_path / "disabled"))
    try:
        assert not run.dir.exists()
        video = Video(data, format="gif")
        assert video.is_valid and video._source is not None
        assert video._source.path is None
        swanlab.log_video(key="clip", data=data, fps=8, format="gif")
        assert not run.dir.exists()
        assert path.read_bytes() == GIF
    finally:
        run.finish()


def test_buffer_assignment_preserves_original_wrapper_source():
    inner = Video(GIF)
    outer = Video(inner)
    replacement = BytesIO(GIF + b"replacement")
    outer.buffer = replacement
    assert outer.buffer is replacement
    assert inner.buffer is not replacement
    assert inner.buffer.getvalue() == GIF
    assert outer.sha256 == hashlib.sha256(GIF + b"replacement").hexdigest()
