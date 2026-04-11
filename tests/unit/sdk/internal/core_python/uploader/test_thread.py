import time
from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.core_python.uploader.thread import Uploader


def test_uploader_defaults_to_class_upload_interval():
    """验证不传 upload_interval 时使用类默认值 UPLOAD_INTERVAL。"""
    with patch("swanlab.sdk.internal.core_python.uploader.thread.timer.Timer") as mock_timer_cls:
        Uploader(auto_start=False)

    _, kwargs = mock_timer_cls.call_args
    assert kwargs["interval"] == Uploader.UPLOAD_INTERVAL


def test_uploader_uses_custom_upload_interval():
    """验证传入自定义 upload_interval 时覆盖类默认值。"""
    with patch("swanlab.sdk.internal.pkg.timer.Timer") as mock_timer_cls:
        Uploader(upload_interval=0.5, auto_start=False)

    _, kwargs = mock_timer_cls.call_args
    assert kwargs["interval"] == 0.5


# ─────────────────── Uploader.start() ───────────────────


def test_uploader_start_reuses_pkg_timer_scheduler():
    """验证 start() 创建 Timer 时参数正确（interval、immediate、name）并启动。"""
    timer = MagicMock()

    with patch("swanlab.sdk.internal.pkg.timer.Timer") as mock_timer_cls:
        mock_timer_cls.return_value = timer

        pool = Uploader(auto_start=False)
        pool.start()

        assert mock_timer_cls.call_count == 1
        _, kwargs = mock_timer_cls.call_args
        assert kwargs["interval"] == Uploader.UPLOAD_INTERVAL
        assert kwargs["immediate"] is True
        assert kwargs["name"] == Uploader.UPLOAD_THREAD_NAME
        timer.start.assert_called_once_with()


def test_uploader_start_is_idempotent():
    """验证 start() 重复调用不会启动第二个 Timer。"""
    timer = MagicMock()

    with patch("swanlab.sdk.internal.pkg.timer.Timer") as mock_timer_cls:
        mock_timer_cls.return_value = timer

        pool = Uploader(auto_start=False)
        pool.start()
        pool.start()

    timer.start.assert_called_once_with()


def test_uploader_start_after_finish_is_noop():
    """验证 finish() 后再次调用 start() 不会重新启动。"""
    pool = Uploader(auto_start=False)
    pool.finish()

    with patch("swanlab.sdk.internal.pkg.timer.Timer") as mock_timer_cls:
        pool.start()

    mock_timer_cls.assert_not_called()


# ─────────────────── Uploader.put() ───────────────────


def test_uploader_put_adds_to_queue(make_config_record):
    """验证 put() 将记录列表放入队列。"""
    pool = Uploader(auto_start=False)
    records = [make_config_record(), make_config_record()]
    pool.put(records)

    assert pool._record_queue.qsize() == 1
    pool.finish()


def test_uploader_put_after_finish_raises(make_config_record):
    """验证 finish() 后调用 put() 抛出 RuntimeError。"""
    pool = Uploader(auto_start=False)
    pool.finish()

    try:
        pool.put([make_config_record()])
        raise AssertionError("put() should raise RuntimeError after finish")
    except RuntimeError as exc:
        assert "already been finished" in str(exc)


# ─────────────────── Uploader._drain_records() ───────────────────


def test_uploader_drain_records_returns_pending_records_once(make_config_record):
    """验证 _drain_records() 取出所有队列记录，再次调用返回空列表。"""
    pool = Uploader(auto_start=False)
    records = [make_config_record(), make_config_record()]
    pool.put(records)

    assert pool._drain_records() == records
    assert pool._drain_records() == []


def test_uploader_drain_records_empty_queue():
    """验证队列为空时 _drain_records() 返回空列表。"""
    pool = Uploader(auto_start=False)

    assert pool._drain_records() == []


# ─────────────────── Uploader.finish() ───────────────────


def test_uploader_finish_is_idempotent(make_config_record):
    """验证 finish() 重复调用只执行一次上传和 Timer 清理。"""
    pool = Uploader(auto_start=False)
    pool.put([make_config_record(), make_config_record()])

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        pool.finish()
        pool.finish()

    assert mock_upload.call_count == 1
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 2
    assert all(record.WhichOneof("record_type") == "config" for record in uploaded)


def test_uploader_finish_cancels_and_joins_timer():
    """验证 finish() 取消并等待 Timer 线程退出。"""
    timer = MagicMock()

    with patch("swanlab.sdk.internal.pkg.timer.Timer") as mock_timer_cls:
        mock_timer_cls.return_value = timer

        pool = Uploader(auto_start=False)
        pool.start()
        pool.finish()

    timer.cancel.assert_called_once_with()
    timer.join.assert_called_once_with(timeout=10)


def test_uploader_finish_drains_and_flushes_remaining(make_scalar_record):
    """验证 finish() 排空队列中的残余记录并通过 flush 上传。"""
    pool = Uploader(auto_start=False)
    pool.put([make_scalar_record(step=1), make_scalar_record(step=2)])

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        pool.finish()

    mock_upload.assert_called_once()
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 2


def test_uploader_finish_passes_upload_callback(make_scalar_record):
    """验证 finish() 通过 Collector 将 upload_callback 透传到上传层。"""
    callback = MagicMock()
    pool = Uploader(upload_interval=0.5, upload_callback=callback, auto_start=False)
    pool.put([make_scalar_record(step=1)])

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        pool.finish()

    mock_upload.assert_called_once()
    # upload_records 第二个参数应为 upload_callback
    assert mock_upload.call_args.kwargs.get("upload_callback") is callback


# ─────────────────── Uploader 集成测试 ───────────────────


def test_uploader_starts_upload_thread_automatically(make_scalar_record):
    """验证 auto_start=True 时自动启动 Timer，put 的记录能在超时后被上传。"""
    with patch.object(Uploader, "UPLOAD_INTERVAL", 0.01):
        pool = Uploader(upload_interval=0.01)

        try:
            with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
                pool.put([make_scalar_record()])

                deadline = time.time() + 1.0
                while mock_upload.call_count == 0 and time.time() < deadline:
                    time.sleep(0.02)

                assert mock_upload.call_count == 1
        finally:
            pool.finish()
