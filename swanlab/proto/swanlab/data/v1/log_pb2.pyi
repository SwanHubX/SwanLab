import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from swanlab.data.v1 import scalar_pb2 as _scalar_pb2
from swanlab.data.v1 import image_pb2 as _image_pb2
from swanlab.data.v1 import audio_pb2 as _audio_pb2
from swanlab.data.v1 import video_pb2 as _video_pb2
from swanlab.data.v1 import text_pb2 as _text_pb2
from swanlab.data.v1 import echarts_pb2 as _echarts_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class LogRecord(_message.Message):
    __slots__ = ("key", "step", "timestamp", "scalar", "images", "audios", "videos", "texts", "echarts")
    KEY_FIELD_NUMBER: _ClassVar[int]
    STEP_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    SCALAR_FIELD_NUMBER: _ClassVar[int]
    IMAGES_FIELD_NUMBER: _ClassVar[int]
    AUDIOS_FIELD_NUMBER: _ClassVar[int]
    VIDEOS_FIELD_NUMBER: _ClassVar[int]
    TEXTS_FIELD_NUMBER: _ClassVar[int]
    ECHARTS_FIELD_NUMBER: _ClassVar[int]
    key: str
    step: int
    timestamp: _timestamp_pb2.Timestamp
    scalar: _scalar_pb2.ScalarValue
    images: _image_pb2.ImageValue
    audios: _audio_pb2.AudioValue
    videos: _video_pb2.VideoValue
    texts: _text_pb2.TextValue
    echarts: _echarts_pb2.EChartsValue
    def __init__(self, key: _Optional[str] = ..., step: _Optional[int] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., scalar: _Optional[_Union[_scalar_pb2.ScalarValue, _Mapping]] = ..., images: _Optional[_Union[_image_pb2.ImageValue, _Mapping]] = ..., audios: _Optional[_Union[_audio_pb2.AudioValue, _Mapping]] = ..., videos: _Optional[_Union[_video_pb2.VideoValue, _Mapping]] = ..., texts: _Optional[_Union[_text_pb2.TextValue, _Mapping]] = ..., echarts: _Optional[_Union[_echarts_pb2.EChartsValue, _Mapping]] = ...) -> None: ...
