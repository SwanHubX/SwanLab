"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 19:17
@description: 音频处理模块
"""

import hashlib
from io import BytesIO
from pathlib import Path
from typing import List, Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.media.audio_pb2 import AudioItem, AudioValue
from swanlab.sdk.internal.context import TransformMediaType
from swanlab.sdk.internal.pkg.fs import safe_write


class Audio(TransformMediaType):
    def __init__(
        self,
        data_or_path: Union["Audio", str, "vendor.np.ndarray"],
        sample_rate: int = 44100,
        caption: Optional[str] = None,
    ):
        """Audio class constructor

        Parameters
        ----------
        data_or_path: str, numpy.ndarray, or Audio
            Path to an audio file, numpy array of audio data (shape: (num_channels, num_frames)),
            or another Audio instance.
        sample_rate: int
            Sample rate of the audio data. Required when input is a numpy array.
        caption: str, optional
            Caption for the audio.
        """
        super().__init__()
        attrs = self._unwrap(data_or_path)
        # 套娃加载
        if attrs:
            self.buffer: BytesIO = attrs["buffer"]
            self.sample_rate: int = attrs["sample_rate"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        if isinstance(data_or_path, str):
            audio_data, sample_rate = vendor.soundfile.read(data_or_path)
            # 转置为 (num_channels, num_frames) 的形式
            audio_data = audio_data.T
        elif isinstance(data_or_path, vendor.np.ndarray):
            sf_support_dtype = [vendor.np.dtype(d) for d in ["float32", "float64", "int16", "int32"]]
            if data_or_path.dtype not in sf_support_dtype:
                raise TypeError(
                    f"Invalid numpy array dtype for audio, supported: {sf_support_dtype}, but got {data_or_path.dtype}"
                )
            # 如果data_or_path是一维, 则reshape为2维
            if len(data_or_path.shape) == 1:
                data_or_path = data_or_path.reshape(1, -1)
            num_channels = data_or_path.shape[0]
            if num_channels not in (1, 2):
                raise TypeError(
                    "Invalid numpy array for audio, support shape is (num_channels, num_frames) with 1 or 2 channels"
                )
            audio_data = data_or_path
        else:
            raise TypeError("Unsupported audio type. Please provide a valid path or numpy array.")

        self.buffer = BytesIO()
        vendor.soundfile.write(self.buffer, audio_data.T, sample_rate, format="wav")
        self.sample_rate = sample_rate
        self.caption = caption

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_AUDIO

    @classmethod
    def build_data_record(cls, *, key: str, step: int, timestamp: Timestamp, data: List[AudioItem]) -> DataRecord:
        return DataRecord(
            key=key, step=step, timestamp=timestamp, type=cls.column_type(), audios=AudioValue(items=data)
        )

    def transform(self, *, step: int, path: Path) -> AudioItem:
        content = self.buffer.getvalue()
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.wav"
        safe_write(path / filename, content, mode="wb")
        return AudioItem(filename=filename, sha256=sha256, size=len(content), caption=self.caption or "")
