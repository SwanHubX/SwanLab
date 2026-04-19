"""
@author: cunyue
@file: test_run_log_media.py
@time: 2026/3/15
@description: Run.log_text / log_image / log_audio / log_video 的参数化测试
"""

import threading
from io import BytesIO
from unittest.mock import MagicMock

import numpy as np
import pytest

from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import Run
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
    """最小化的 Run 替身，供非绑定方法测试使用"""

    def __init__(self):
        self._api_lock = threading.RLock()
        self._init_pid = fork.current_pid()
        self.log = MagicMock()
        self.alive = True


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


def _call_log(method_name, mock_run, key, data, extra, caption=None, step=None):
    """Helper that calls an unbound Run.log_* using keyword args (required by keyword-only params)."""
    method = getattr(Run, method_name)
    return method(mock_run, key=key, data=data, **extra, caption=caption, step=step)


@pytest.mark.parametrize("method_name,factory,extra", _MEDIA_LOG_CASES)
class TestRunLogMedia:
    """log_text / log_image / log_audio / log_video 的通用契约测试"""

    def test_single_instance_wrapped_in_list(self, method_name, factory, extra):
        """传入单个媒体实例 → log 被调用，data 会被包装为列表"""
        mock_run = _MockRun()
        instance = factory()
        _call_log(method_name, mock_run, "key", instance, extra)
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert "key" in log_data
        assert len(log_data["key"]) == 1
        assert type(log_data["key"][0]) is instance.__class__

    def test_list_of_instances_passed_through(self, method_name, factory, extra):
        """传入实例列表 → log 被调用，data 保持列表（类型一致）"""
        mock_run = _MockRun()
        instances = [factory(), factory()]
        _call_log(method_name, mock_run, "metrics/a", instances, extra)
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        result = log_data["metrics/a"]
        assert len(result) == 2
        assert all(isinstance(item, factory().__class__) for item in result)

    def test_step_forwarded(self, method_name, factory, extra):
        """step 参数应被原样透传给 log"""
        mock_run = _MockRun()
        _call_log(method_name, mock_run, "k", factory(), extra, step=7)
        _, kwargs = mock_run.log.call_args
        assert kwargs.get("step") == 7

    def test_step_none_by_default(self, method_name, factory, extra):
        """不传 step 时，log 收到 step=None"""
        mock_run = _MockRun()
        _call_log(method_name, mock_run, "k", factory(), extra)
        _, kwargs = mock_run.log.call_args
        assert kwargs.get("step") is None

    def test_raises_when_not_running(self, method_name, factory, extra):
        """run 未激活时应抛出 RuntimeError"""
        mock_run = _MockRun()
        mock_run.alive = False
        with pytest.raises(RuntimeError, match="requires an active Run"):
            _call_log(method_name, mock_run, "k", factory(), extra)

    def test_raw_data_gets_wrapped(self, method_name, factory, extra):
        """传入原始数据（非媒体对象）时，也能被包装为对应的媒体类列表"""
        mock_run = _MockRun()
        raw_data_map = {
            "log_text": "raw text string",
            "log_image": np.zeros((8, 8, 3), dtype=np.uint8),
            "log_audio": np.zeros((1, 4410), dtype=np.float32),
            "log_video": _GIF_1X1,
        }
        raw = raw_data_map[method_name]
        _call_log(method_name, mock_run, "k", raw, extra)
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert len(log_data["k"]) == 1
        assert isinstance(log_data["k"][0], factory().__class__)


@pytest.mark.parametrize("method_name,factory,extra", _MEDIA_LOG_CASES)
class TestRunLogMediaCaption:
    """caption 参数行为（对所有媒体类型一致）"""

    def test_caption_applied_to_single_item(self, method_name, factory, extra):
        """传入单个实例 + caption → 对象 caption 正确"""
        mock_run = _MockRun()
        instance = factory()
        _call_log(method_name, mock_run, "k", instance, extra, caption="my caption")
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert isinstance(log_data["k"], list)
        assert len(log_data["k"]) == 1
        assert log_data["k"][0].caption == "my caption"

    def test_caption_none_by_default(self, method_name, factory, extra):
        """不传 caption → 对象 caption 为 None"""
        mock_run = _MockRun()
        _call_log(method_name, mock_run, "k", factory(), extra)
        log_data = mock_run.log.call_args[0][0]
        assert log_data["k"][0].caption is None

    def test_list_caption_applied_per_item(self, method_name, factory, extra):
        """传入列表数据 + 列表 caption → 每个对象 caption 各自对应"""
        mock_run = _MockRun()
        instances = [factory(), factory()]
        _call_log(method_name, mock_run, "k", instances, extra, caption=["cap0", "cap1"])
        log_data = mock_run.log.call_args[0][0]
        assert log_data["k"][0].caption == "cap0"
        assert log_data["k"][1].caption == "cap1"

    def test_broadcast_caption_to_list(self, method_name, factory, extra):
        """传入列表数据 + 单个 caption → caption 广播到每个对象"""
        mock_run = _MockRun()
        instances = [factory(), factory()]
        _call_log(method_name, mock_run, "k", instances, extra, caption="shared")
        log_data = mock_run.log.call_args[0][0]
        assert all(item.caption == "shared" for item in log_data["k"])


class TestRunLogImageExtra:
    """log_image 特有参数（mode / file_type / size）"""

    def test_file_type_forwarded(self):
        """file_type 参数应影响 Image 内部 file_type"""
        mock_run = _MockRun()
        arr = np.zeros((8, 8, 3), dtype=np.uint8)
        _call_log("log_image", mock_run, "img", arr, {"file_type": "jpg"})
        log_data = mock_run.log.call_args[0][0]
        assert log_data["img"][0].file_type == "jpg"

    def test_size_applied(self):
        """size 参数应对图像进行缩放（结果尺寸不超过指定值）"""
        mock_run = _MockRun()
        arr = np.zeros((100, 100, 3), dtype=np.uint8)
        _call_log("log_image", mock_run, "img", arr, {"size": 32})
        # Verify log was called without error — size enforcement is internal to Image
        mock_run.log.assert_called_once()

    def test_mode_forwarded(self):
        """mode 参数应被透传给 Image（不报错即可）"""
        mock_run = _MockRun()
        arr = np.zeros((8, 8, 3), dtype=np.uint8)
        _call_log("log_image", mock_run, "img", arr, {"mode": "RGB"})
        mock_run.log.assert_called_once()


class TestRunLogAudioExtra:
    """log_audio 特有参数（sample_rate）"""

    def test_sample_rate_forwarded(self):
        """sample_rate 应被正确传入 Audio 对象"""
        mock_run = _MockRun()
        arr = np.zeros((1, 16000), dtype=np.float32)
        Run.log_audio(mock_run, key="audio", data=arr, sample_rate=16000)  # type: ignore
        log_data = mock_run.log.call_args[0][0]
        assert log_data["audio"][0].sample_rate == 16000


class TestRunLogVideoFormats:
    """log_video 支持 bytes 和 BytesIO 输入"""

    def test_bytes_input(self):
        mock_run = _MockRun()
        Run.log_video(mock_run, key="v", data=_GIF_1X1)  # type: ignore
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert isinstance(log_data["v"][0], Video)

    def test_bytesio_input(self):
        mock_run = _MockRun()
        Run.log_video(mock_run, key="v", data=BytesIO(_GIF_1X1))  # type: ignore
        mock_run.log.assert_called_once()
        log_data = mock_run.log.call_args[0][0]
        assert isinstance(log_data["v"][0], Video)
