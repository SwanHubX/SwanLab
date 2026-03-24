"""
@author: CaddiesNew
@file: save_manager.py
@time: 2026/3/24 14:10
@description: swanlab.save() 的上传/监听管理器
"""

import hashlib
import math
import os
import shutil
import threading
import time
from dataclasses import dataclass
from io import BytesIO
from queue import Queue
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

import requests
from requests.exceptions import RequestException

from swanlab.core_python import get_client
from swanlab.core_python.api.file_service import (
    MULTIPART_THRESHOLD,
    PART_SIZE,
    complete_multipart,
    complete_upload,
    prepare_multipart,
    prepare_upload,
    upload_file,
)
from swanlab.log import swanlog

from .utils import SaveFile


@dataclass
class _WatchEntry:
    file: SaveFile
    signature: Optional[Tuple[int, int]]


def _file_signature(path: str) -> Optional[Tuple[int, int]]:
    try:
        stat = os.stat(path)
    except FileNotFoundError:
        return None
    return stat.st_mtime_ns, stat.st_size


def _file_md5(path: str, chunk_size: int = 1024 * 1024) -> str:
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if chunk == b"":
                break
            md5.update(chunk)
    return md5.hexdigest()


def _remove_path(path: str) -> None:
    if not os.path.lexists(path):
        return
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.unlink(path)


def _iter_files(files: Union[SaveFile, Iterable[SaveFile]]) -> Iterable[SaveFile]:
    if isinstance(files, SaveFile):
        return [files]
    return files


def _extract_upload_url(payload: Dict[str, object]) -> str:
    for key in ("uploadUrl", "upload_url", "url", "presignedUrl", "presigned_url"):
        value = payload.get(key)
        if isinstance(value, str) and value != "":
            return value
    raise ValueError("Upload URL is missing in prepare response.")


def _extract_upload_id(payload: Dict[str, object]) -> Optional[str]:
    for key in ("uploadId", "upload_id"):
        value = payload.get(key)
        if isinstance(value, str) and value != "":
            return value
    return None


def _extract_part_size(payload: Dict[str, object]) -> int:
    for key in ("partSize", "part_size"):
        value = payload.get(key)
        if isinstance(value, int) and value > 0:
            return value
    return PART_SIZE


def _extract_part_urls(payload: Dict[str, object]) -> List[Tuple[int, str]]:
    parts = payload.get("parts")
    if isinstance(parts, list):
        resolved = []
        for index, part in enumerate(parts, start=1):
            if not isinstance(part, dict):
                raise ValueError("Multipart prepare response contains invalid part data.")
            number = part.get("partNumber", part.get("part_number", index))
            resolved.append((int(number), _extract_upload_url(part)))
        return sorted(resolved, key=lambda item: item[0])

    urls = payload.get("uploadUrls", payload.get("upload_urls", payload.get("urls")))
    if isinstance(urls, list):
        resolved = []
        for index, url in enumerate(urls, start=1):
            if not isinstance(url, str) or url == "":
                raise ValueError("Multipart prepare response contains invalid upload URL.")
            resolved.append((index, url))
        return resolved

    raise ValueError("Multipart upload URLs are missing in prepare response.")


def _upload_part(url: str, buffer: BytesIO, max_retries: int = 3) -> Optional[str]:
    with requests.Session() as session:
        for attempt in range(1, max_retries + 1):
            try:
                buffer.seek(0)
                response = session.put(
                    url,
                    data=buffer,
                    headers={"Content-Type": "application/octet-stream"},
                    timeout=30,
                )
                response.raise_for_status()
                etag = response.headers.get("ETag") or response.headers.get("etag")
                return etag.strip('"') if isinstance(etag, str) else None
            except RequestException:
                swanlog.warning(f"Upload attempt {attempt} failed for multipart URL: {url}")
                if attempt == max_retries:
                    raise
                time.sleep(2 ** (attempt - 1))


class FileUploadManager:
    """
    每个 Run 独立的文件处理线程。
    cloud 模式下执行上传；local/offline 模式下复制到 run/files；disabled 模式下跳过。
    """

    _STOP = object()

    def __init__(self, mode: str, file_dir: str):
        self._mode = mode
        self._file_dir = file_dir
        self._client = get_client() if mode == "cloud" else None
        self._exp_id: Optional[str] = None
        self._queue = Queue()
        self._closed = False
        self._thread = threading.Thread(target=self._worker, name="swanlab-save-upload", daemon=True)
        self._thread.start()

    def submit(self, files: Union[SaveFile, Iterable[SaveFile]]) -> None:
        if self._closed:
            return
        for file in _iter_files(files):
            self._queue.put(file)

    def join(self) -> None:
        self._queue.join()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._queue.put(self._STOP)
        self._thread.join(timeout=5)

    def _worker(self) -> None:
        while True:
            item = self._queue.get()
            try:
                if item is self._STOP:
                    return
                self._do_upload(item)
            except Exception as e:
                swanlog.warning(f"Failed to save file {item.name}: {e}")
            finally:
                self._queue.task_done()

    def _do_upload(self, file: SaveFile) -> None:
        if not os.path.exists(file.source_path):
            swanlog.warning(f"Save file not found, skip it: {file.source_path}")
            return
        if self._mode == "disabled":
            return
        if self._mode in ("local", "offline"):
            self._copy_to_run_dir(file)
            return
        if self._mode != "cloud":
            raise ValueError(f"Unsupported save mode: {self._mode}")
        self._upload_to_cloud(file)

    def _copy_to_run_dir(self, file: SaveFile) -> None:
        dst_path = file.target_path
        if os.path.abspath(file.source_path) == os.path.abspath(dst_path):
            return
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        _remove_path(dst_path)
        shutil.copy2(file.source_path, dst_path)

    def _upload_to_cloud(self, file: SaveFile) -> None:
        assert self._client is not None, "Client must be initialized before uploading save files."
        exp_id = self._get_exp_id()
        size = os.path.getsize(file.source_path)
        if size < MULTIPART_THRESHOLD:
            self._upload_single_part(file, size, exp_id)
        else:
            self._upload_multipart(file, size, exp_id)

    def _get_exp_id(self) -> str:
        if self._exp_id is not None:
            return self._exp_id
        assert self._client is not None, "Client must be initialized before uploading save files."
        try:
            self._exp_id = self._client.exp_id
        except (AssertionError, AttributeError) as e:
            raise AssertionError("Experiment id must be initialized before uploading save files.") from e
        return self._exp_id

    def _upload_single_part(self, file: SaveFile, size: int, exp_id: str) -> None:
        prepared = prepare_upload(
            self._client,
            exp_id,
            [{"name": file.name, "size": size, "md5": _file_md5(file.source_path)}],
        )
        if len(prepared) == 0:
            raise ValueError("Prepare upload returned an empty result.")
        upload_url = _extract_upload_url(prepared[0])
        with open(file.source_path, "rb") as f:
            upload_file(url=upload_url, buffer=BytesIO(f.read()))
        complete_upload(self._client, exp_id, [file.name])

    def _upload_multipart(self, file: SaveFile, size: int, exp_id: str) -> None:
        part_count = max(1, math.ceil(size / PART_SIZE))
        prepared = prepare_multipart(self._client, exp_id, file.name, size, part_count)
        upload_id = _extract_upload_id(prepared)
        part_size = _extract_part_size(prepared)
        upload_urls = _extract_part_urls(prepared)
        parts = []

        with open(file.source_path, "rb") as f:
            for part_number, upload_url in upload_urls:
                chunk = f.read(part_size)
                if chunk == b"":
                    break
                etag = _upload_part(upload_url, BytesIO(chunk))
                part_info: Dict[str, object] = {"partNumber": part_number}
                if etag is not None:
                    part_info["etag"] = etag
                parts.append(part_info)

        complete_multipart(self._client, exp_id, file.name, parts, upload_id=upload_id)


class DirWatcher:
    """
    仅用于 live policy 的目录监听器。
    cloud 模式下会在 run/files 中创建 symlink 便于保留文件结构；其余模式仅负责监听源文件变化。
    """

    def __init__(self, mode: str, file_dir: str, on_change: Callable[[SaveFile], None], interval: float = 1.0):
        self._mode = mode
        self._file_dir = file_dir
        self._on_change = on_change
        self._interval = interval
        self._entries: Dict[str, _WatchEntry] = {}
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def watch(self, files: Union[SaveFile, Iterable[SaveFile]]) -> None:
        with self._lock:
            for file in _iter_files(files):
                if self._mode == "cloud":
                    self._create_live_link(file)
                self._entries[file.name] = _WatchEntry(file=file, signature=_file_signature(file.source_path))
        self._ensure_started()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop_event.set()
        self._thread.join(timeout=max(1.0, self._interval + 1.0))

    def _ensure_started(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._thread = threading.Thread(target=self._poll, name="swanlab-save-watch", daemon=True)
        self._thread.start()

    def _poll(self) -> None:
        while not self._stop_event.wait(self._interval):
            with self._lock:
                entries = list(self._entries.items())
            for key, entry in entries:
                signature = _file_signature(entry.file.source_path)
                if signature == entry.signature:
                    continue
                with self._lock:
                    latest = self._entries.get(key)
                    if latest is None:
                        continue
                    latest.signature = signature
                    current_file = latest.file
                if signature is None:
                    continue
                try:
                    self._on_change(current_file)
                except Exception as e:
                    swanlog.warning(f"Failed to handle live save change for {current_file.name}: {e}")

    def _create_live_link(self, file: SaveFile) -> None:
        link_path = file.target_path
        if os.path.abspath(link_path) == os.path.abspath(file.source_path):
            return
        os.makedirs(os.path.dirname(link_path), exist_ok=True)
        if os.path.lexists(link_path):
            if os.path.islink(link_path):
                current_target = os.readlink(link_path)
                if not os.path.isabs(current_target):
                    current_target = os.path.abspath(os.path.join(os.path.dirname(link_path), current_target))
                if os.path.abspath(current_target) == os.path.abspath(file.source_path):
                    return
            _remove_path(link_path)
        try:
            target = os.path.relpath(file.source_path, os.path.dirname(link_path))
            os.symlink(target, link_path)
        except (NotImplementedError, OSError):
            swanlog.debug(f"Failed to create symlink for live save: {file.source_path}")


__all__ = [
    "SaveFile",
    "FileUploadManager",
    "DirWatcher",
]
