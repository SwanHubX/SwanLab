"""
@author: cunyue
@file: v0.py
@time: 2025/6/20 13:57
@description: Json 序列化协议上传器，此为使用 protobuf 的前身
此版本无法通过 core 上传数据
"""

from swanlab.core_python import (
    ColumnModel,
    MediaModel,
    FileModel,
    LogModel,
    ScalarModel,
)
from swanlab.proto.v0 import Log, Column, Media, Runtime, Scalar
from .py import PythonTransfer


class ProtoV0Transfer(PythonTransfer):

    def transfer_scalar(self, data: Scalar) -> ScalarModel:
        return data.to_scalar_model()

    def transfer_column(self, data: Column) -> ColumnModel:
        return data.to_column_model()

    def transfer_media(self, data: Media) -> MediaModel:
        return data.to_media_model(self.media_dir)

    def transfer_file(self, data: Runtime) -> FileModel:
        return data.to_file_model(self.file_dir)

    def transfer_log(self, data: Log) -> LogModel:
        return data.to_log_model()

    def publish_column(self, data: Column):
        super().publish_column(data)

    def publish_media(self, data: Media):
        super().publish_media(data)

    def publish_scalar(self, data: Scalar):
        super().publish_scalar(data)

    def publish_file(self, data: Runtime):
        super().publish_file(data)

    def publish_log(self, data: Log):
        super().publish_log(data)
