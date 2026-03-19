"""
@author: caddiesnew
@file: model.py
@time: 2026/3/19
@description: 上传模型定义，基于 protobuf Record 的适配层
"""

from datetime import datetime
from enum import Enum
from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


class UploadType(Enum):
    """上传类型枚举，对应 Record 的 oneof record_type 分类"""

    SCALAR_METRIC = "scalar"
    MEDIA_METRIC = "media"
    COLUMN = "column"
    LOG = "log"
    FILE = "file"


# 媒体类型的 value oneof 名称集合
_MEDIA_VALUE_FIELDS = {"images", "audios", "videos", "texts", "echarts"}
# 文件类型的 record_type 名称集合
_FILE_RECORD_TYPES = {"metadata", "requirements", "conda", "config"}


def classify_record(record: Record) -> Optional[UploadType]:
    """
    根据 Record 的 oneof record_type 判断上传类型。
    返回 None 表示该 Record 不需要走上传线程（如 run/finish 生命周期消息）。
    """
    record_type = record.WhichOneof("record_type")
    if record_type is None:
        return None
    if record_type == "metric":
        value_type = record.metric.WhichOneof("value")
        if value_type == "scalar":
            return UploadType.SCALAR_METRIC
        if value_type in _MEDIA_VALUE_FIELDS:
            return UploadType.MEDIA_METRIC
        return None
    if record_type == "column":
        return UploadType.COLUMN
    if record_type == "console":
        return UploadType.LOG
    if record_type in _FILE_RECORD_TYPES:
        return UploadType.FILE
    # run / finish 等生命周期消息不走上传线程
    return None


class FileModel:
    """
    运行时文件信息上传模型，聚合 metadata/requirements/conda/config 四种 Record。
    """

    def __init__(
        self,
        requirements: Optional[str] = None,
        metadata: Optional[dict] = None,
        config: Optional[dict] = None,
        conda: Optional[str] = None,
    ):
        self.requirements = requirements
        self.metadata = metadata
        self.config = config
        self.conda = conda
        self.create_time = datetime.now()

    @classmethod
    def create(cls, file_models: List["FileModel"]) -> "FileModel":
        """
        比较若干个 FileModel，获取最新的 FileModel，并且保证其内部属性不为 None。
        """
        file_models = sorted(file_models, key=lambda x: x.create_time, reverse=True)
        lr, lm, lc, lo = None, None, None, None
        for fm in file_models:
            lr = fm.requirements if lr is None else lr
            lm = fm.metadata if lm is None else lm
            lc = fm.config if lc is None else lc
            lo = fm.conda if lo is None else lo
            if lr is not None and lm is not None and lc is not None and lo is not None:
                break
        return FileModel(lr, lm, lc, lo)

    def to_dict(self) -> dict:
        """序列化，删除为 None 的字段。"""
        d = {
            "requirements": self.requirements,
            "metadata": self.metadata,
            "config": self.config,
            "conda": self.conda,
        }
        return {k: v for k, v in d.items() if v is not None}

    @property
    def empty(self) -> bool:
        """是否为空。"""
        return all(getattr(self, attr) is None for attr in ("requirements", "metadata", "config", "conda"))


__all__ = [
    "UploadType",
    "FileModel",
    "classify_record",
]
