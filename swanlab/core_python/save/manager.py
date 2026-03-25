"""
@author: CaddiesNew
@file: manager.py
@time: 2026/3/24 14:10
@description: swanlab.save() 的上传/监听管理器
"""

import hashlib
import math
import mimetypes
import os
import shutil
import threading
from concurrent.futures import Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from io import BytesIO
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

from swanlab.core_python.api.experiment import (
    MULTIPART_THRESHOLD,
    PART_SIZE,
    complete_multipart,
    complete_upload,
    prepare_multipart,
    prepare_upload,
)
from swanlab.core_python.api.experiment.utils import (
    extract_part_size,
    extract_part_urls,
    extract_upload_id,
    extract_upload_url,
)
from swanlab.core_python.api.service import upload_file
from swanlab.core_python.client import get_client
from swanlab.log import swanlog

from .model import SaveFile


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


def _guess_mime_type(path: str) -> Optional[str]:
    mime_type, _ = mimetypes.guess_type(path)
    return mime_type


def _build_prepare_file_payload(file: SaveFile, size: int) -> Dict[str, object]:
    payload: Dict[str, object] = {
        "name": file.name,
        "size": size,
        "md5": _file_md5(file.source_path),
    }
    mime_type = _guess_mime_type(file.source_path)
    if mime_type is not None:
        payload["mimeType"] = mime_type
    return payload


def _upload_buffers(buffers: List[Tuple[str, BytesIO]]) -> None:
    if len(buffers) == 0:
        return

    failed_buffers: List[Tuple[str, BytesIO]] = []
    last_error: Optional[Exception] = None
    max_workers = min(10, len(buffers))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures: List[Tuple[Future, str, BytesIO]] = []
        for url, buffer in buffers:
            try:
                future = executor.submit(upload_file, url=url, buffer=buffer)
                futures.append((future, url, buffer))
            except RuntimeError:
                failed_buffers.append((url, buffer))

        for future, url, buffer in futures:
            try:
                future.result()
            except Exception as e:
                swanlog.warning(f"Failed to upload {url}: {e}, will retry...")
                failed_buffers.append((url, buffer))

    if len(failed_buffers):
        swanlog.debug(f"Retrying failed save buffers: {len(failed_buffers)}")
        for url, buffer in failed_buffers:
            try:
                upload_file(url=url, buffer=buffer)
            except Exception as e:
                swanlog.error(f"Failed to upload {url}: {e}")
                last_error = e

    if last_error is not None:
        raise last_error


class FileUploadManager:
    """
    每个 Run 独立的文件处理线程。
    cloud 模式下执行上传；local/offline 模式下复制到 run/files；disabled 模式下跳过。
    """

    def __init__(self, mode: str, file_dir: str, max_workers: int = 4):
        self._mode = mode
        self._file_dir = file_dir
        self._client = get_client() if mode == "cloud" else None
        self._exp_id: Optional[str] = None
        self._closed = False
        self._lock = threading.Lock()
        self._executor = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="swanlab-save-upload")
        self._futures: Dict[Future, str] = {}

    def submit(self, files: Union[SaveFile, Iterable[SaveFile]]) -> None:
        if self._closed:
            return
        for file in _iter_files(files):
            try:
                future = self._executor.submit(self._do_upload, file)
            except RuntimeError:
                if self._closed:
                    return
                try:
                    self._do_upload(file)
                except Exception as e:
                    swanlog.warning(f"Failed to save file {file.name}: {e}")
                continue

            with self._lock:
                self._futures[future] = file.name
            future.add_done_callback(self._consume_future)

    def join(self) -> None:
        while True:
            with self._lock:
                futures = list(self._futures.keys())
            if len(futures) == 0:
                return
            wait(futures)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self.join()
        self._executor.shutdown(wait=True)

    def _consume_future(self, future: Future) -> None:
        with self._lock:
            file_name = self._futures.pop(future, "<unknown>")
        try:
            future.result()
        except Exception as e:
            swanlog.warning(f"Failed to save file {file_name}: {e}")

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
        file_payload = _build_prepare_file_payload(file, size)
        prepared = prepare_upload(
            self._client,
            exp_id,
            [file_payload],
        )
        if len(prepared) == 0:
            raise ValueError("Prepare upload returned an empty result.")
        upload_url = extract_upload_url(prepared[0])
        with open(file.source_path, "rb") as f:
            upload_file(url=upload_url, buffer=f)
        complete_upload(self._client, exp_id, [file.name])

    def _upload_multipart(self, file: SaveFile, size: int, exp_id: str) -> None:
        part_count = max(1, math.ceil(size / PART_SIZE))
        file_payload = _build_prepare_file_payload(file, size)
        prepared = prepare_multipart(
            self._client,
            exp_id,
            file.name,
            size,
            part_count,
            md5=str(file_payload["md5"]),
            mime_type=file_payload.get("mimeType"),
        )
        upload_id = extract_upload_id(prepared)
        if upload_id is None:
            raise ValueError("Upload ID is missing in multipart prepare response.")
        part_size = extract_part_size(prepared, PART_SIZE)
        upload_urls = extract_part_urls(prepared)
        buffers: List[Tuple[str, BytesIO]] = []

        with open(file.source_path, "rb") as f:
            for _part_number, upload_url in upload_urls:
                chunk = f.read(part_size)
                if chunk == b"":
                    break
                buffers.append((upload_url, BytesIO(chunk)))

        _upload_buffers(buffers)
        complete_multipart(self._client, exp_id, file.name, upload_id=upload_id)


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
