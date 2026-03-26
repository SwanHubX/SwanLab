from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple

FileSignature = Tuple[int, int]


class SaveFileState(str, Enum):
    UPLOADED = "UPLOADED"
    FAILED = "FAILED"


@dataclass
class WatchSaveFileModel:
    source_path: str
    name: str
    target_path: str
    signature: Optional[FileSignature] = None

    def prepare_request(
        self,
        *,
        size: int,
        md5: str,
        mime_type: Optional[str] = None,
        count: Optional[int] = None,
    ) -> Dict[str, object]:
        payload: Dict[str, object] = {
            "path": self.name,
            "size": size,
            "md5": md5,
        }
        if mime_type is not None:
            payload["mimeType"] = mime_type
        if count is not None:
            payload["count"] = count
        return payload

    def complete_request(
        self,
        *,
        state: SaveFileState = SaveFileState.UPLOADED,
    ) -> Dict[str, object]:
        return {
            "path": self.name,
            "state": state.value,
        }

    def complete_multipart_request(
        self,
        *,
        upload_id: str,
        parts: List[Dict[str, object]],
    ) -> Dict[str, object]:
        return {
            "path": self.name,
            "uploadId": upload_id,
            "parts": parts,
        }


__all__ = [
    "FileSignature",
    "SaveFileState",
    "WatchSaveFileModel",
]
