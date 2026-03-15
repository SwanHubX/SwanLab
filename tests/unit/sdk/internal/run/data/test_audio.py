"""
@author: cunyue
@file: test_audio.py
@time: 2026/3/15
@description: 测试音频转换模块 Audio
"""

import hashlib

import numpy as np
import pytest
import soundfile as sf
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.media.audio_pb2 import AudioItem
from swanlab.sdk.internal.run.transforms.audio import Audio

# ---------------------------------- Fixtures ----------------------------------


@pytest.fixture
def mono_array():
    """单声道 numpy 数组，shape=(1, 4410)"""
    return np.zeros((1, 4410), dtype=np.float32)


@pytest.fixture
def stereo_array():
    """双声道 numpy 数组，shape=(2, 4410)"""
    return np.zeros((2, 4410), dtype=np.float32)


@pytest.fixture
def wav_file(tmp_path):
    """生成一个临时 WAV 文件路径"""
    path = tmp_path / "test.wav"
    data = np.zeros((4410, 1), dtype=np.float32)  # soundfile 写入为 (frames, channels)
    sf.write(str(path), data, 44100)
    return str(path)


# ---------------------------------- 构造测试 ----------------------------------


class TestAudioInit:
    def test_from_mono_array(self, mono_array):
        """从单声道 numpy 数组构造"""
        audio = Audio(mono_array, sample_rate=44100)
        assert audio.sample_rate == 44100
        assert audio.caption is None
        assert len(audio.buffer.getvalue()) > 0

    def test_from_stereo_array(self, stereo_array):
        """从双声道 numpy 数组构造"""
        audio = Audio(stereo_array, sample_rate=22050)
        assert audio.sample_rate == 22050
        assert len(audio.buffer.getvalue()) > 0

    def test_from_1d_array_auto_reshape(self):
        """一维数组自动 reshape 为单声道"""
        arr = np.zeros(4410, dtype=np.float32)
        audio = Audio(arr, sample_rate=44100)
        assert len(audio.buffer.getvalue()) > 0

    def test_from_wav_path(self, wav_file):
        """从 WAV 文件路径构造"""
        audio = Audio(wav_file)
        assert len(audio.buffer.getvalue()) > 0

    def test_caption_stored(self, mono_array):
        """caption 被正确保存"""
        audio = Audio(mono_array, caption="test caption")
        assert audio.caption == "test caption"

    def test_caption_none_by_default(self, mono_array):
        """caption 默认为 None"""
        audio = Audio(mono_array)
        assert audio.caption is None


class TestAudioInitErrors:
    def test_invalid_dtype_raises(self):
        """不支持的 dtype 应抛出 TypeError"""
        arr = np.zeros((1, 4410), dtype=np.int8)
        with pytest.raises(TypeError, match="Invalid numpy array dtype"):
            Audio(arr)

    def test_invalid_channels_raises(self):
        """超过2声道应抛出 TypeError"""
        arr = np.zeros((3, 4410), dtype=np.float32)
        with pytest.raises(TypeError, match="1 or 2 channels"):
            Audio(arr)

    def test_unsupported_type_raises(self):
        """传入不支持的类型应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported audio type"):
            Audio(12345)  # type: ignore

    def test_invalid_path_raises(self):
        """无效路径应抛出 ValueError"""
        with pytest.raises(Exception):
            Audio("/nonexistent/path/audio.wav")


# ---------------------------------- 套娃加载测试 ----------------------------------


class TestAudioNesting:
    def test_wrap_audio_copies_buffer(self, mono_array):
        """套娃加载复用内层 buffer"""
        inner = Audio(mono_array, sample_rate=44100, caption="inner")
        outer = Audio(inner)
        assert outer.buffer is inner.buffer
        assert outer.sample_rate == inner.sample_rate
        assert outer.caption == "inner"

    def test_outer_caption_overrides_inner(self, mono_array):
        """外层 caption 优先级高于内层"""
        inner = Audio(mono_array, caption="inner")
        outer = Audio(inner, caption="outer")
        assert outer.caption == "outer"

    def test_inner_caption_used_when_outer_none(self, mono_array):
        """外层 caption 为 None 时使用内层 caption"""
        inner = Audio(mono_array, caption="inner")
        outer = Audio(inner, caption=None)
        assert outer.caption == "inner"


# ---------------------------------- column_type 测试 ----------------------------------


class TestAudioColumnType:
    def test_column_type(self):
        assert Audio.column_type() == ColumnType.COLUMN_TYPE_AUDIO


# ---------------------------------- build_data_record 测试 ----------------------------------


class TestAudioBuildDataRecord:
    def test_build_data_record_structure(self, mono_array, tmp_path):
        """build_data_record 返回正确结构的 DataRecord"""
        audio = Audio(mono_array)
        item = audio.transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Audio.build_data_record(key="loss", step=1, timestamp=ts, data=[item])

        assert record.key == "loss"
        assert record.step == 1
        assert record.type == ColumnType.COLUMN_TYPE_AUDIO
        assert len(record.audios.items) == 1
        assert record.audios.items[0].filename == item.filename
        assert record.audios.items[0].sha256 == item.sha256

    def test_build_data_record_multiple_items(self, mono_array, stereo_array, tmp_path):
        """build_data_record 支持多个 AudioItem"""
        a1 = Audio(mono_array).transform(step=1, path=tmp_path)
        a2 = Audio(stereo_array).transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Audio.build_data_record(key="k", step=1, timestamp=ts, data=[a1, a2])
        assert len(record.audios.items) == 2


# ---------------------------------- transform 测试 ----------------------------------


class TestAudioTransform:
    def test_transform_returns_audio_item(self, mono_array, tmp_path):
        """transform 返回 AudioItem"""
        item = Audio(mono_array).transform(step=1, path=tmp_path)
        assert isinstance(item, AudioItem)

    def test_transform_sha256_correct(self, mono_array, tmp_path):
        """AudioItem.sha256 与落盘文件内容一致"""
        item = Audio(mono_array).transform(step=1, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        assert item.sha256 == hashlib.sha256(content).hexdigest()

    def test_transform_size_correct(self, mono_array, tmp_path):
        """AudioItem.size 与落盘文件字节数一致"""
        item = Audio(mono_array).transform(step=1, path=tmp_path)
        assert item.size == len((tmp_path / item.filename).read_bytes())

    def test_transform_caption_empty_when_none(self, mono_array, tmp_path):
        """caption 为 None 时，AudioItem.caption 为空字符串"""
        item = Audio(mono_array).transform(step=1, path=tmp_path)
        assert item.caption == ""

    def test_transform_caption_preserved(self, mono_array, tmp_path):
        """caption 有值时，AudioItem.caption 正确"""
        item = Audio(mono_array, caption="hello").transform(step=1, path=tmp_path)
        assert item.caption == "hello"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
