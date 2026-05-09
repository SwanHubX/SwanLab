"""
@author: caddiesnew
@file: save.py
@time: 2026/4/30
@description: 文件保存相关API类型定义
"""

from typing import List, Literal, TypedDict


class SaveFileEntry(TypedDict):
    """Save 文件入口"""

    path: str
    # 源路径
    source_path: str
    size: int
    md5: str
    mime_type: str


class PrepareSaveFile(TypedDict):
    """prepare 请求中的单个文件描述"""

    path: str
    size: int
    md5: str
    mimeType: str


class PrepareSaveResponse(TypedDict):
    """prepare 接口响应，预签名 URL 列表，顺序与请求中 files 顺序一致"""

    urls: List[str]


class CompleteSaveFile(TypedDict):
    """complete 请求中的单个文件"""

    path: str
    state: Literal["UPLOADED", "FAILED"]


class PrepareMultipartSaveFile(TypedDict):
    """prepare-multipart 请求中的单个文件描述"""

    path: str
    size: int
    md5: str
    mimeType: str
    count: int


class MultipartPart(TypedDict):
    """分片上传的单个分片信息"""

    partNumber: int
    url: str


class MultipartFilePayload(TypedDict):
    """prepare-multipart 返回的单个文件的分片上传信息"""

    uploadId: str
    parts: List[MultipartPart]


class PrepareMultipartSaveResponse(TypedDict):
    """prepare-multipart 接口响应，顺序与请求中 files 顺序一致"""

    files: List[MultipartFilePayload]


class CompletedPart(TypedDict):
    """分片上传完成时报告的单个分片状态"""

    partNumber: int
    etag: str


class CompletedMultipartSaveFile(TypedDict):
    """complete-multipart 请求中的单个文件信息"""

    path: str
    parts: List[CompletedPart]
    uploadId: str
