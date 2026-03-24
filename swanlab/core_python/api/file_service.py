"""
@author: CaddiesNew
@file: file_service.py
@time: 2026/3/24 14:10
@description: save() 所需的文件上传 API
"""

from typing import TYPE_CHECKING, Dict, List, Optional

from .service import upload_file

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


MULTIPART_THRESHOLD = 100 * 1024 * 1024
PART_SIZE = 5 * 1024 * 1024


def _unwrap_payload(data):
    if isinstance(data, dict) and "data" in data and isinstance(data["data"], (dict, list)):
        return data["data"]
    return data


def prepare_upload(client: "Client", exp_id: str, files: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """
    创建普通文件上传任务，返回预签名上传地址列表。
    """
    if len(files) == 0:
        return []
    data, _ = client.post(f"/experiment/{exp_id}/files/prepare", files)
    payload = _unwrap_payload(data)
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        if any(isinstance(payload.get(key), str) for key in ("uploadUrl", "upload_url", "url", "presignedUrl")):
            return [payload]
        for key in ("files", "items", "list"):
            value = payload.get(key)
            if isinstance(value, list):
                return value
    raise ValueError("Unexpected response from file prepare API.")


def complete_upload(client: "Client", exp_id: str, names: List[str]) -> None:
    """
    标记普通文件上传完成。
    """
    if len(names) == 0:
        return
    client.post(f"/experiment/{exp_id}/files/complete", [{"name": name} for name in names])


def prepare_multipart(
    client: "Client",
    exp_id: str,
    name: str,
    size: int,
    part_count: int,
) -> Dict[str, object]:
    """
    创建分片上传任务，返回上传地址和上传上下文。
    """
    data, _ = client.post(
        f"/experiment/{exp_id}/files/prepare-multipart",
        {"name": name, "size": size, "partCount": part_count},
    )
    payload = _unwrap_payload(data)
    if not isinstance(payload, dict):
        raise ValueError("Unexpected response from multipart prepare API.")
    return payload


def complete_multipart(
    client: "Client",
    exp_id: str,
    name: str,
    parts: List[Dict[str, object]],
    upload_id: Optional[str] = None,
) -> None:
    """
    标记分片上传完成，并通知后端执行合并。
    """
    payload = {"name": name, "parts": parts}
    if upload_id is not None:
        payload["uploadId"] = upload_id
    client.post(f"/experiment/{exp_id}/files/complete-multipart", payload)


__all__ = [
    "MULTIPART_THRESHOLD",
    "PART_SIZE",
    "prepare_upload",
    "complete_upload",
    "prepare_multipart",
    "complete_multipart",
    "upload_file",
]
