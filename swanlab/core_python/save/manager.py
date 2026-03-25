"""swanlab.save() 的上传/监听管理器"""

import math
import os
import shutil
import threading
from concurrent.futures import Future, ThreadPoolExecutor, wait
from io import BytesIO
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

from swanlab.core_python.api.experiment import (
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

from .model import WatchSaveFileModel
from .utils import compute_md5, file_signature, guess_mime_type

# 分片阈值 / 大小
MULTIPART_THRESHOLD: int = 100 * 1024 * 1024
PART_SIZE = 10 * 1024 * 1024

def _remove_path(path: str) -> None:
    if not os.path.lexists(path):
        return
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.unlink(path)


def _iter_files(files: Union[WatchSaveFileModel, Iterable[WatchSaveFileModel]]) -> Iterable[WatchSaveFileModel]:
    return [files] if isinstance(files, WatchSaveFileModel) else files


def _upload_buffers(buffers: List[Tuple[str, BytesIO]]) -> None:
    """并发上传多个分片"""
    if not buffers:
        return

    failed = []
    with ThreadPoolExecutor(max_workers=min(10, len(buffers))) as executor:
        futures = [(executor.submit(upload_file, url=url, buffer=buf), url, buf) for url, buf in buffers]
        for future, url, buf in futures:
            try:
                future.result()
            except Exception as e:
                swanlog.warning(f"Upload failed: {e}, will retry")
                failed.append((url, buf))

    for url, buf in failed:
        upload_file(url=url, buffer=buf)


class FileUploadManager:
    """文件上传管理器，支持 cloud/local/disabled 模式"""

    def __init__(self, mode: str, file_dir: str, max_workers: int = 4):
        self._mode = mode
        self._file_dir = file_dir
        self._client = get_client() if mode == "cloud" else None
        self._exp_id: Optional[str] = None
        self._closed = False
        self._lock = threading.Lock()
        self._executor = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="swanlab-save")
        self._futures: Dict[Future, str] = {}

    def submit(self, files: Union[WatchSaveFileModel, Iterable[WatchSaveFileModel]]) -> None:
        if self._closed:
            return
        for file in _iter_files(files):
            try:
                future = self._executor.submit(self._do_upload, file)
            except RuntimeError:
                if not self._closed:
                    self._do_upload(file)
                continue

            with self._lock:
                self._futures[future] = file.name
            future.add_done_callback(self._on_done)

    def join(self) -> None:
        while True:
            with self._lock:
                futures = list(self._futures.keys())
            if not futures:
                return
            wait(futures)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self.join()
        self._executor.shutdown(wait=True)

    def _on_done(self, future: Future) -> None:
        with self._lock:
            name = self._futures.pop(future, "<unknown>")
        try:
            future.result()
        except Exception as e:
            swanlog.warning(f"Save failed for {name}: {e}")

    def _do_upload(self, file: WatchSaveFileModel) -> None:
        if not os.path.exists(file.source_path):
            swanlog.warning(f"File not found: {file.source_path}")
            return
        if self._mode == "disabled":
            return
        if self._mode in ("local", "offline"):
            self._copy_local(file)
            return
        if self._mode == "cloud":
            self._upload_cloud(file)

    def _copy_local(self, file: WatchSaveFileModel) -> None:
        if os.path.abspath(file.source_path) == os.path.abspath(file.target_path):
            return
        os.makedirs(os.path.dirname(file.target_path), exist_ok=True)
        _remove_path(file.target_path)
        shutil.copy2(file.source_path, file.target_path)

    def _upload_cloud(self, file: WatchSaveFileModel) -> None:
        exp_id = self._get_exp_id()
        size = os.path.getsize(file.source_path)
        if size < MULTIPART_THRESHOLD:
            self._upload_single(file, size, exp_id)
        else:
            self._upload_multipart(file, size, exp_id)

    def _get_exp_id(self) -> str:
        if self._exp_id:
            return self._exp_id
        assert self._client is not None
        exp_id = self._client.exp_id
        assert exp_id is not None
        self._exp_id = exp_id
        return self._exp_id

    def _upload_single(self, file: WatchSaveFileModel, size: int, exp_id: str) -> None:
        assert self._client is not None
        payload = file.prepare_payload(
            size=size,
            md5=compute_md5(file.source_path),
            mime_type=guess_mime_type(file.source_path),
        )
        prepared = prepare_upload(self._client, exp_id, [payload])
        upload_url = extract_upload_url(prepared[0])
        with open(file.source_path, "rb") as f:
            upload_file(url=upload_url, buffer=BytesIO(f.read()))
        complete_upload(self._client, exp_id, [file.complete_payload()])

    def _upload_multipart(self, file: WatchSaveFileModel, size: int, exp_id: str) -> None:
        assert self._client is not None
        part_count = max(1, math.ceil(size / PART_SIZE))
        payload = file.prepare_payload(
            size=size,
            md5=compute_md5(file.source_path),
            mime_type=guess_mime_type(file.source_path),
            count=part_count,
        )
        prepared = prepare_multipart(self._client, exp_id, payload)
        upload_id = extract_upload_id(prepared)
        assert upload_id is not None
        part_size = extract_part_size(prepared, PART_SIZE)
        upload_urls = extract_part_urls(prepared)

        buffers = []
        with open(file.source_path, "rb") as f:
            for _, url in upload_urls:
                chunk = f.read(part_size)
                if chunk != b"":
                    buffers.append((url, BytesIO(chunk)))

        _upload_buffers(buffers)
        complete_multipart(self._client, exp_id, file.complete_payload(upload_id=upload_id))


class DirWatcher:
    """目录监听器，用于 live policy"""

    def __init__(self, mode: str, file_dir: str, on_change: Callable[[WatchSaveFileModel], None], interval: float = 1.0):
        self._mode = mode
        self._file_dir = file_dir
        self._on_change = on_change
        self._interval = interval
        self._entries: Dict[str, WatchSaveFileModel] = {}
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def watch(self, files: Union[WatchSaveFileModel, Iterable[WatchSaveFileModel]]) -> None:
        with self._lock:
            for file in _iter_files(files):
                if self._mode == "cloud":
                    self._create_symlink(file)
                file.signature = file_signature(file.source_path)
                self._entries[file.name] = file
        self._start()

    def stop(self) -> None:
        if self._thread:
            self._stop_event.set()
            self._thread.join(timeout=self._interval + 1.0)

    def _start(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._thread = threading.Thread(target=self._poll, name="swanlab-save-watch", daemon=True)
        self._thread.start()

    def _poll(self) -> None:
        while not self._stop_event.wait(self._interval):
            with self._lock:
                entries = list(self._entries.items())
            for name, file in entries:
                sig = file_signature(file.source_path)
                if sig == file.signature or sig is None:
                    continue
                with self._lock:
                    if name not in self._entries:
                        continue
                    self._entries[name].signature = sig
                    current = self._entries[name]
                try:
                    self._on_change(current)
                except Exception as e:
                    swanlog.warning(f"Live save change failed for {name}: {e}")

    def _create_symlink(self, file: WatchSaveFileModel) -> None:
        if os.path.abspath(file.target_path) == os.path.abspath(file.source_path):
            return
        os.makedirs(os.path.dirname(file.target_path), exist_ok=True)
        if os.path.islink(file.target_path):
            target = os.readlink(file.target_path)
            if not os.path.isabs(target):
                target = os.path.abspath(os.path.join(os.path.dirname(file.target_path), target))
            if os.path.abspath(target) == os.path.abspath(file.source_path):
                return
        _remove_path(file.target_path)
        try:
            rel_target = os.path.relpath(file.source_path, os.path.dirname(file.target_path))
            os.symlink(rel_target, file.target_path)
        except (NotImplementedError, OSError):
            pass


__all__ = ["FileUploadManager", "DirWatcher"]
