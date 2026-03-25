from dataclasses import dataclass
from typing import Literal, Optional, Tuple
from enum import Enum

@dataclass
class SaveFileModel:
    """文件上传模型，对应后端 prepare 接口的 file 对象"""
    name: str
    size: int
    md5: str
    count: Optional[int]
    mimeType: Optional[str] = None

    def to_dict(self, multipart: bool = False):
        d = {"name": self.name, "size": self.size, "md5": self.md5}
        if self.mimeType:
            d["mimeType"] = self.mimeType
        if multipart:
            d["count"] = self.count if self.count is not None else 1
        return d


class SaveFileState(Enum):
    UPLOADED = "UPLOADED"
    FAILED = "FAILED"



@dataclass
class SaveFileStateModel:
    """文件状态模型，对应后端 complete 接口的 file 对象"""
    name: str
    state: SaveFileState

    def to_dict(self):
        return {"name": self.name, "state": self.state.value}


@dataclass
class WatchSaveFileModel:
    """监听文件模型，包含本地路径和上传信息"""
    fid: int
    source_path: str
    name: str
    target_path: str
    signature: Optional[tuple] = None

    def to_save_file_model(self, size: int, md5: str, mime_type: Optional[str] = None, count: Optional[int] = None) -> SaveFileModel:
        return SaveFileModel(name=self.name, size=size, md5=md5, mimeType=mime_type, count=count)


__all__ = ["SaveFileModel", "SaveFileState", "SaveFileStateModel", "WatchSaveFileModel"]
