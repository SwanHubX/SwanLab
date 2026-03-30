import threading
import time
from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.core_python.uploader.collector import Collector

# ─────────────────── upload() ───────────────────


def test_upload_delegates_to_upload_records(make_scalar_record):
    """验证 upload() 将 pending 列表传递给 upload_records。"""
    collector = Collector()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.upload(records)

    mock_upload.assert_called_once_with(records, upload_callback=None)


def test_upload_passes_callback(make_scalar_record):
    """验证 upload() 将构造时传入的 upload_callback 透传给 upload_records。"""
    callback = MagicMock()
    collector = Collector(upload_callback=callback)
    records = [make_scalar_record(step=1)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.upload(records)

    mock_upload.assert_called_once_with(records, upload_callback=callback)


def test_upload_skips_empty_list():
    """验证 upload() 传入空列表时不调用 upload_records。"""
    collector = Collector()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.upload([])

    mock_upload.assert_not_called()


# ─────────────────── submit() ───────────────────


def test_submit_accumulates_and_uploads_when_interval_elapsed(make_scalar_record):
    """验证 submit() 在超过 upload_interval 后聚合并上传累积的记录。"""
    collector = Collector(upload_interval=0.0)
    records = [make_scalar_record(step=i) for i in range(3)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.submit(records)

    mock_upload.assert_called_once()
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 3


def test_submit_uploads_immediately_when_zero_interval_and_clock_does_not_advance(make_scalar_record):
    """验证 upload_interval=0.0 时会立即上传。"""
    records = [make_scalar_record(step=1)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.time.time", side_effect=[100.0, 100.0]):
        collector = Collector(upload_interval=0.0)

        with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
            collector.submit(records)

    mock_upload.assert_called_once()
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 1


def test_submit_accumulates_without_upload_when_interval_not_elapsed(make_scalar_record):
    """验证 submit() 在未超过 upload_interval 时只聚合不上传。"""
    collector = Collector(upload_interval=999.0)
    records = [make_scalar_record(step=1)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.submit(records)

    mock_upload.assert_not_called()
    assert len(collector.container) == 1


def test_submit_empty_records_is_noop():
    """验证 submit() 传入空列表时无任何副作用。"""
    collector = Collector()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.submit([])

    mock_upload.assert_not_called()
    assert collector.container == []


def test_submit_accumulates_across_calls_then_uploads(make_scalar_record):
    """验证多次 submit() 调用累积记录，直到超过 interval 才一次性上传。"""
    collector = Collector(upload_interval=0.0)
    records_1 = [make_scalar_record(step=1)]
    records_2 = [make_scalar_record(step=2)]

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.submit(records_1)
        collector.submit(records_2)

    assert mock_upload.call_count == 2
    # 第二次调用时 container 里只有 records_2，因为第一次已清空
    second_uploaded = mock_upload.call_args_list[1].args[0]
    assert len(second_uploaded) == 1


def test_submit_rolls_back_on_upload_error(make_scalar_record):
    """验证 submit() 上传失败时，pending 记录被回滚到 container 头部。"""
    collector = Collector(upload_interval=0.0)
    # 先放入一条记录模拟已积累
    collector.container.append(make_scalar_record(step=0))
    new_records = [make_scalar_record(step=1)]

    with patch(
        "swanlab.sdk.internal.core_python.uploader.collector.upload_records", side_effect=RuntimeError("network down")
    ), patch("swanlab.sdk.internal.core_python.uploader.collector.console") as mock_console:
        collector.submit(new_records)

    mock_console.error.assert_called_once()
    assert "network down" in mock_console.error.call_args[0][0]
    # 回滚后 container 应为 pending + 新增 = [old, new]
    assert len(collector.container) == 2
    assert collector.container[0].metric.step == 0
    assert collector.container[1].metric.step == 1


# ─────────────────── flush() ───────────────────


def test_flush_flushes_all_accumulated_records(make_scalar_record):
    """验证 flush() 无条件将所有累积记录一次性上传。"""
    collector = Collector(upload_interval=999.0)
    collector.container.append(make_scalar_record(step=1))

    new_records = [make_scalar_record(step=2)]
    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.flush(new_records)

    mock_upload.assert_called_once()
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 2
    assert collector.container == []


def test_flush_empty_list_no_upload():
    """验证 flush() 传入空列表且 container 也为空时不触发上传。"""
    collector = Collector()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.flush([])

    mock_upload.assert_not_called()


def test_flush_rolls_back_on_upload_error(make_scalar_record):
    """验证 flush() 上传失败时，pending 记录被回滚到 container。"""
    collector = Collector()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]

    with patch(
        "swanlab.sdk.internal.core_python.uploader.collector.upload_records", side_effect=ValueError("bad data")
    ), patch("swanlab.sdk.internal.core_python.uploader.collector.console") as mock_console:
        collector.flush(records)

    mock_console.error.assert_called_once()
    assert "bad data" in mock_console.error.call_args[0][0]
    assert len(collector.container) == 2


def test_flush_waits_for_lock_release(make_scalar_record):
    """验证 flush() 在锁被持有时阻塞等待，释放后继续执行。"""
    collector = Collector()
    record = make_scalar_record(step=1)
    done = threading.Event()

    collector._lock.acquire()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        worker = threading.Thread(target=lambda: (collector.flush([record]), done.set()))
        worker.start()

        time.sleep(0.05)
        assert not done.is_set()
        mock_upload.assert_not_called()

        collector._lock.release()
        worker.join(timeout=1)

    assert done.is_set()
