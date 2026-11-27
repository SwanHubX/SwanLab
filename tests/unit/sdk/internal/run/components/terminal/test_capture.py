"""StreamCapture：透传优先，console.skip_capture() 期间不触发采集回调（告警自循环防线）。"""

import io
import os
import sys
import threading

import pytest

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run.components.terminal.capture import StreamCapture


def _install(monkeypatch, name):
    stream = io.StringIO()
    monkeypatch.setattr(sys, name, stream)
    captured = []
    capture = StreamCapture(name, captured.append, os.getpid())
    capture.install()
    return capture, stream, captured


def test_write_is_passed_through_and_captured(monkeypatch):
    capture, stream, captured = _install(monkeypatch, "stdout")
    try:
        print("hello")
    finally:
        capture.uninstall()

    assert stream.getvalue() == "hello\n"
    # print() 分两次 write：正文 + 换行
    assert captured == ["hello", "\n"]


def test_stderr_is_captured_independently(monkeypatch):
    capture, stream, captured = _install(monkeypatch, "stderr")
    try:
        print("boom", file=sys.stderr)
    finally:
        capture.uninstall()

    assert stream.getvalue() == "boom\n"
    assert captured == ["boom", "\n"]


def test_console_skip_capture_keeps_output_but_skips_capture(monkeypatch):
    """console.skip_capture()：告警照常上屏，但不作为用户日志上报，且异常后自动复位。"""
    capture, stream, captured = _install(monkeypatch, "stdout")
    try:
        console.info("Tracking run with swanlab version 0.11.0-dev")
        with console.skip_capture():
            console.warning("Dropped 3 event(s): consumer is not running")
        with pytest.raises(RuntimeError, match="boom"):
            with console.skip_capture():
                with console.skip_capture():
                    assert console.capture_skipped() is True
                assert console.capture_skipped() is True
                raise RuntimeError("boom")
        assert console.capture_skipped() is False
        console.info("Resumed normal capture")
    finally:
        capture.uninstall()

    screen = " ".join(stream.getvalue().split())
    replayed = " ".join("".join(captured).split())
    assert "Dropped 3 event(s)" in screen
    assert "Tracking run with swanlab version 0.11.0-dev" in replayed
    assert "Resumed normal capture" in replayed
    assert "Dropped 3 event(s)" not in replayed


def test_skip_capture_does_not_leak_to_other_threads():
    """skip_capture 仅作用于当前上下文（线程），不能误伤其他线程的正常采集。"""
    seen = []

    with console.skip_capture():
        assert console.capture_skipped() is True
        thread = threading.Thread(target=lambda: seen.append(console.capture_skipped()))
        thread.start()
        thread.join()

    assert seen == [False]
    assert console.capture_skipped() is False
