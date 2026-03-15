"""
@author: cunyue
@file: test_run_log_media.py
@time: 2026/3/15
@description: SwanLabRun.log_text / log_image / log_audio / log_video 的参数化测试
"""

import threading
from unittest.mock import MagicMock

import numpy as np
import pytest

from swanlab.sdk.internal.run import SwanLabRun
from swanlab.sdk.internal.run.transforms.audio import Audio
from swanlab.sdk.internal.run.transforms.image import Image
from swanlab.sdk.internal.run.transforms.text import Text
from swanlab.sdk.internal.run.transforms.video import Video

# 最小合法 GIF89a（1×1 像素）
_GIF_1X1 = (
    b"GIF89a\x01\x00\x01\x00\x80\x00\x00\xff\xff\xff\x00\x00\x00"
    b"!\xf9\x04\x00\x00\x00\x00\x00,\x00\x00\x00\x00\x01\x00\x01\x00\x00\x02\x02D\x01\x00;"
)


def _make_text() -> Text:
    return Text(content="hello world")


def _make_audio() -> Audio:
    return Audio(np.zeros((1, 4410), dtype=np.float32), sample_rate=44100)


def _make_image() -> Image:
    return Image(np.zeros((10, 10, 3), dtype=np.uint8))


def _make_video() -> Video:
    return Video(_GIF_1X1)


# ──────────────────────────────────────────────
# Mock run
# ──────────────────────────────────────────────


class _MockRun:
    """最小化的 SwanLabRun 替身，供非绑定方法测试使用"""

    def __init__(self):
        self._state = "running"
        self._api_lock = threading.RLock()
        self.log = MagicMock()


# ──────────────────────────────────────────────
# 参数化用例
# ──────────────────────────────────────────────

# (method_name, media_factory, extra_kwargs_for_method)
_MEDIA_LOG_CASES = [
    pytest.param("log_text", _make_text, {}, id="log_text"),
    pytest.param("log_image", _make_image, {}, id="log_image"),
    pytest.param("log_audio", _make_audio, {"sample_rate": 44100}, id="log_audio"),
    pytest.param("log_video", _make_video, {}, id="log_video"),
]


def _call_step(method_name, mock_run, key, data, extra, caption=None, step=None):
    method = getattr(SwanLabRun, method_name)
    return method(mock_run, key, data, **extra, caption=caption, step=step)


@pytest.mark.parametrize("method_name,factory,extra", _MEDIA_LOG_CASES)
class TestRunLogMedia:
    """log_text / log_image / log_audio / log_video 的通用契约测试
    没有log_echarts和log_object3D，因为他们的操作比较复杂
    """

    def test_single_instance_wrapped_in_list(self, method_name, factory, extra):
        """传入单个媒体实例 → log 被调用，data 会被包装"""
        mock_run = _MockRun()
        instance = factory()
        _call_step(method_name, mock_run, "key", instance, extra)
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert "key" in log_data
        assert len(log_data["key"]) == 1
        assert log_data["key"] != [instance]
        assert type(log_data["key"][0]) is instance.__class__

    def test_list_of_instances_passed_through(self, method_name, factory, extra):
        """传入实例列表 → log 被调用，data 保持列表不变"""
        mock_run = _MockRun()
        instances = [factory(), factory()]
        _call_step(method_name, mock_run, "metrics/a", instances, extra)
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert log_data["metrics/a"] != instances
        assert all(isinstance(item, factory().__class__) for item in log_data["metrics/a"])

    def test_step_forwarded(self, method_name, factory, extra):
        """step 参数应被原样透传给 log"""
        mock_run = _MockRun()
        _call_step(method_name, mock_run, "k", factory(), extra, step=7)
        _, kwargs = mock_run.log.call_args
        assert kwargs.get("step") == 7

    def test_step_none_by_default(self, method_name, factory, extra):
        """不传 step 时，log 收到 step=None"""
        mock_run = _MockRun()
        _call_step(method_name, mock_run, "k", factory(), extra)
        _, kwargs = mock_run.log.call_args
        assert kwargs.get("step") is None

    def test_raises_when_not_running(self, method_name, factory, extra):
        """run 未激活时应抛出 RuntimeError"""
        mock_run = _MockRun()
        mock_run._state = "finished"
        with pytest.raises(RuntimeError, match="requires an active SwanLabRun"):
            _call_step(method_name, mock_run, "k", factory(), extra)


def _call_caption(method_name, mock_run, key, data, extra, caption=None, step=None):
    method = getattr(SwanLabRun, method_name)
    return method(mock_run, key, data, **extra, caption=caption, step=step)


@pytest.mark.parametrize("method_name,factory,extra", _MEDIA_LOG_CASES)
class TestRunLogMediaCaption:
    """caption 参数行为（对所有媒体类型一致）"""

    def test_caption_applied_to_raw_data(self, method_name, factory, extra):
        """传入原始数据 + caption → normalize 后的对象 caption 正确"""
        mock_run = _MockRun()
        instance = factory()
        _call_caption(method_name, mock_run, "k", instance, extra, caption="my caption")
        # We can't inspect the caption easily without calling transform, so just
        # verify log was called with a list
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert isinstance(log_data["k"], list)
        assert len(log_data["k"]) == 1
        assert log_data["k"][0].caption == "my caption"
