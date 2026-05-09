"""
@author: caddiesnew
@file: save.py
@time: 2026/4/30
@description: 文件保存相关API：prepare / complete / prepare-multipart / complete-multipart
"""

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.save import (
    CompletedMultipartSaveFile,
    CompleteSaveFile,
    PrepareMultipartSaveFile,
    PrepareMultipartSaveResponse,
    PrepareSaveFile,
    PrepareSaveResponse,
)


def prepare_save_files(experiment_id: str, *, files: list[PrepareSaveFile]) -> PrepareSaveResponse:
    """批量准备单文件上传，返回预签名 URL 列表。

    响应 data 结构：{"urls": ["https://...", ...]}
    """
    resp = client.post(f"/experiment/{experiment_id}/files/prepare", {"files": files})
    return resp.data


def complete_save_files(experiment_id: str, *, files: list[CompleteSaveFile]) -> None:
    """批量完成单文件上传，报告上传状态。

    状态码：201 成功，404 文件不存在。
    """
    client.post(f"/experiment/{experiment_id}/files/complete", {"files": files})


def prepare_multipart_save(
    experiment_id: str, *, files: list[PrepareMultipartSaveFile]
) -> PrepareMultipartSaveResponse:
    """批量准备分片上传，返回分片预签名 URL。

    响应 data 结构：{"files": [{"uploadId": "...", "parts": [{"partNumber": 1, "url": "..."}]}]}
    """
    resp = client.post(f"/experiment/{experiment_id}/files/prepare-multipart", {"files": files})
    return resp.data


def complete_multipart_save(experiment_id: str, *, files: list[CompletedMultipartSaveFile]) -> None:
    """批量完成分片上传，报告分片 ETag。

    状态码：201 成功，400 参数错误，404 文件不存在。
    """
    client.post(f"/experiment/{experiment_id}/files/complete-multipart", {"files": files})
