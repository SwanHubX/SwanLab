import time
from unittest.mock import patch

from swanlab.sdk.internal.core_python.transport.thread import UploadWarningThrottle


def test_warn_throttles_within_interval():
    """间隔内重复 warn() 不重复打印，超过间隔后再次打印。

    INTERVAL 设为 1.0s 而非极小值，因为 Windows 上 time.monotonic() 分辨率约 15ms，
    过小的间隔会导致连续两次 warn() 被判定为"已超过间隔"，节流失效。
    """
    t = UploadWarningThrottle()
    t.INTERVAL = 1.0
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning") as mock:
        t.warn()
        t.warn()
    mock.assert_called_once()

    time.sleep(1.1)
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning") as mock:
        t.warn()
    mock.assert_called_once()


def test_reset_clears_throttle_and_prints_recovery():
    """失败状态下 reset() 打印恢复消息，之后 warn() 不再被节流。"""
    t = UploadWarningThrottle()
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning"):
        t.warn()

    with patch("swanlab.sdk.internal.core_python.transport.thread.console.info") as mock_info:
        t.reset()
    mock_info.assert_called_once()

    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning") as mock_warn:
        t.warn()
    mock_warn.assert_called_once()


def test_reset_without_failure_no_recovery():
    """未进入失败状态时 reset() 不打印恢复消息。"""
    t = UploadWarningThrottle()
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.info") as mock_info:
        t.reset()
    mock_info.assert_not_called()
