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
    extract_part_urls,
    extract_upload_id,
)
from swanlab.core_python.api.service import upload_file
from swanlab.core_python.client import get_client
from swanlab.core_python.utils.timer import Timer
from swanlab.log import swanlog

from .model import SaveFileModel, SaveFileState
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


def _iter_files(
    files: Union[SaveFileModel, Iterable[SaveFileModel]],
) -> Iterable[SaveFileModel]:
    return [files] if isinstance(files, SaveFileModel) else files


def _normalize_etag(etag: Optional[str]) -> str:
    if not isinstance(etag, str) or etag == "":
        raise ValueError("Multipart upload response is missing ETag.")
    return etag.strip('"')


def _upload_buffers(buffers: List[Tuple[int, str, BytesIO]]) -> List[Dict[str, object]]:
    """并发上传多个分片并返回完成合并所需的 etag 信息"""
    if not buffers:
        return []

    failed: List[Tuple[int, str, BytesIO]] = []
    completed: Dict[int, Dict[str, object]] = {}
    with ThreadPoolExecutor(max_workers=min(10, len(buffers))) as executor:
        futures = [
            (executor.submit(upload_file, url=url, buffer=buf), part_number, url, buf)
            for part_number, url, buf in buffers
        ]
        for future, part_number, url, buf in futures:
            try:
                completed[part_number] = {
                    "partNumber": part_number,
                    "etag": _normalize_etag(future.result()),
                }
            except Exception as e:
                swanlog.warning(f"Upload failed: {e}, will retry")
                failed.append((part_number, url, buf))

    for part_number, url, buf in failed:
        etag = upload_file(url=url, buffer=buf)
        completed[part_number] = {
            "partNumber": part_number,
            "etag": _normalize_etag(etag),
        }

    return [completed[part_number] for part_number in sorted(completed)]


class FileUploadManager:
    """文件上传管理器，支持 cloud/local/disabled 模式"""

    def __init__(self, mode: str, file_dir: str, max_workers: int = 4):
        self._mode = mode
        self._file_dir = file_dir
        self._client = get_client() if mode == "cloud" else None
        self._exp_id: Optional[str] = None
        self._closed = False
        self._lock = threading.Lock()
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="swanlab-save"
        )
        self._futures: Dict[Future, str] = {}

    def submit(self, files: Union[SaveFileModel, Iterable[SaveFileModel]]) -> None:
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

    def _do_upload(self, file: SaveFileModel) -> None:
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

    def _copy_local(self, file: SaveFileModel) -> None:
        if os.path.abspath(file.source_path) == os.path.abspath(file.target_path):
            return
        os.makedirs(os.path.dirname(file.target_path), exist_ok=True)
        _remove_path(file.target_path)
        shutil.copy2(file.source_path, file.target_path)

    def _upload_cloud(self, file: SaveFileModel) -> None:
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

    def _upload_single(self, file: SaveFileModel, size: int, exp_id: str) -> None:
        assert self._client is not None
        mime_type = guess_mime_type(file.source_path)
        payload = file.prepare_request(
            size=size,
            md5=compute_md5(file.source_path),
            mime_type=mime_type,
        )
        urls = prepare_upload(self._client, exp_id, [payload])
        try:
            with open(file.source_path, "rb") as f:
                upload_file(url=urls[0], buffer=BytesIO(f.read()), mime_type=mime_type)
            complete_upload(self._client, exp_id, [file.complete_request()])
        except Exception as e:
            complete_upload(
                self._client,
                exp_id,
                [file.complete_request(state=SaveFileState.FAILED)],
            )

    def _upload_multipart(self, file: SaveFileModel, size: int, exp_id: str) -> None:
        assert self._client is not None
        part_count = max(1, math.ceil(size / PART_SIZE))
        mime_type = guess_mime_type(file.source_path)
        payload = file.prepare_request(
            size=size,
            md5=compute_md5(file.source_path),
            mime_type=mime_type,
            count=part_count,
        )
        prepared = prepare_multipart(self._client, exp_id, payload)
        upload_id = extract_upload_id(prepared)
        assert upload_id is not None
        upload_urls = extract_part_urls(prepared)

        buffers = []
        with open(file.source_path, "rb") as f:
            for part_number, url in upload_urls:
                chunk = f.read(PART_SIZE)
                if chunk:
                    buffers.append((part_number, url, BytesIO(chunk)))

        parts = _upload_buffers(buffers)
        complete_multipart(
            self._client,
            exp_id,
            file.complete_multipart_request(upload_id=upload_id, parts=parts),
        )


class DirWatcher:
    """目录监听器，用于 live policy"""

    def __init__(
        self,
        mode: str,
        file_dir: str,
        on_change: Callable[[SaveFileModel], None],
        interval: float = 1.0,
    ):
        self._mode = mode
        self._file_dir = file_dir
        self._on_change = on_change
        self._interval = interval
        self._entries: Dict[str, SaveFileModel] = {}
        self._lock = threading.Lock()
        self._timer = Timer(task=self._poll_task, interval=self._interval)

    def watch(self, files: Union[SaveFileModel, Iterable[SaveFileModel]]) -> None:
        with self._lock:
            for file in _iter_files(files):
                if self._mode == "cloud":
                    self._create_symlink(file)
                file.signature = file_signature(file.source_path)
                self._entries[file.name] = file
        self._start()

    def stop(self) -> None:
        self._timer.cancel()
        self._timer.join(timeout=self._interval + 1.0)

    def _start(self) -> None:
        self._timer.run()

    def _poll_task(self) -> None:
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

    def _create_symlink(self, file: SaveFileModel) -> None:
        if os.path.abspath(file.target_path) == os.path.abspath(file.source_path):
            return
        os.makedirs(os.path.dirname(file.target_path), exist_ok=True)
        if os.path.islink(file.target_path):
            target = os.readlink(file.target_path)
            if not os.path.isabs(target):
                target = os.path.abspath(
                    os.path.join(os.path.dirname(file.target_path), target)
                )
            if os.path.abspath(target) == os.path.abspath(file.source_path):
                return
        _remove_path(file.target_path)
        try:
            rel_target = os.path.relpath(
                file.source_path, os.path.dirname(file.target_path)
            )
            os.symlink(rel_target, file.target_path)
        except (NotImplementedError, OSError):
            pass


__all__ = ["FileUploadManager", "DirWatcher"]
