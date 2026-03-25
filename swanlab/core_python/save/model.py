from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple

FileSignature = Tuple[int, int]


class SaveFileState(str, Enum):
    UPLOADED = "UPLOADED"
    FAILED = "FAILED"


@dataclass(frozen=True, slots=True)
class SaveFilePayload:
    name: str
    size: Optional[int] = None
    md5: Optional[str] = None
    mime_type: Optional[str] = None
    count: Optional[int] = None
    upload_id: Optional[str] = None
    state: Optional[SaveFileState] = None

    def to_dict(self) -> dict:
        payload = {"name": self.name}
        if self.size is not None:
            payload["size"] = self.size
        if self.md5 is not None:
            payload["md5"] = self.md5
        if self.mime_type is not None:
            payload["mimeType"] = self.mime_type
        if self.count is not None:
            payload["count"] = self.count
        if self.upload_id is not None:
            payload["uploadId"] = self.upload_id
        if self.state is not None:
            payload["state"] = self.state.value
        return payload


@dataclass(slots=True)
class WatchSaveFileModel:
    source_path: str
    name: str
    target_path: str
    signature: Optional[FileSignature] = None

    def prepare_payload(
        self,
        *,
        size: int,
        md5: str,
        mime_type: Optional[str] = None,
        count: Optional[int] = None,
    ) -> SaveFilePayload:
        return SaveFilePayload(
            name=self.name,
            size=size,
            md5=md5,
            mime_type=mime_type,
            count=count,
        )

    def complete_payload(
        self,
        *,
        upload_id: Optional[str] = None,
        state: SaveFileState = SaveFileState.UPLOADED,
    ) -> SaveFilePayload:
        return SaveFilePayload(
            name=self.name,
            upload_id=upload_id,
            state=state,
        )


__all__ = [
    "FileSignature",
    "SaveFileState",
    "SaveFilePayload",
    "WatchSaveFileModel",
]
