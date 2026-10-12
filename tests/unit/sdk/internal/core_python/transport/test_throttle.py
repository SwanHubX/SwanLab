import time
from unittest.mock import patch

from swanlab.sdk.internal.core_python.transport.worker import UploadWarningThrottle


def test_warn_throttles_within_interval():
    """间隔内调用 warn() 抑制重复告警，超过间隔后输出告警。"""
    t = UploadWarningThrottle()
    t.INTERVAL = 1.0
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.warning") as mock:
        t.warn()
        t.warn()
    mock.assert_called_once()

    time.sleep(1.1)
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.warning") as mock:
        t.warn()
    mock.assert_called_once()


def test_reset_clears_throttle_and_prints_recovery():
    """失败状态下调用 reset() 输出恢复信息并重置节流状态。"""
    t = UploadWarningThrottle()
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.warning"):
        t.warn()

    with patch("swanlab.sdk.internal.core_python.transport.worker.console.info") as mock_info:
        t.reset()
    mock_info.assert_called_once()

    with patch("swanlab.sdk.internal.core_python.transport.worker.console.warning") as mock_warn:
        t.warn()
    mock_warn.assert_called_once()


def test_reset_without_failure_no_recovery():
    """未处于失败状态时 reset() 不输出恢复信息。"""
    t = UploadWarningThrottle()
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.info") as mock_info:
        t.reset()
    mock_info.assert_not_called()
