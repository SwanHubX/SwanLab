"""swanlab.save() 的上传/监听管理器"""

import math
import os
import pathlib
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
from .progress import _SaveProgress
from .utils import compute_md5, copy_file, file_signature, guess_mime_type, same_path

# 分片阈值 / 大小
MULTIPART_THRESHOLD: int = 100 * 1024 * 1024
PART_SIZE = 10 * 1024 * 1024

# 单文件上传大小上限 (50 GB)
MAX_FILE_SIZE: int = 50 * 1024 * 1024 * 1024


def _iter_files(
    files: Union[SaveFileModel, Iterable[SaveFileModel]],
) -> List[SaveFileModel]:
    return [files] if isinstance(files, SaveFileModel) else list(files)


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
        futures = []
        for part_number, url, buf in buffers:
            try:
                future = executor.submit(upload_file, url=url, buffer=buf)
                futures.append((future, part_number, url, buf))
            except RuntimeError:
                failed.append((part_number, url, buf))
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
        self._progress: Optional[_SaveProgress] = None

    def submit(self, files: Union[SaveFileModel, Iterable[SaveFileModel]]) -> None:
        if self._closed:
            return
        file_list = list(_iter_files(files))
        # 仅在 cloud 模式下展示上传进度
        if self._mode == "cloud" and file_list:
            if self._progress is None:
                self._progress = _SaveProgress(len(file_list))
            else:
                self._progress.add(len(file_list))
        submitted: Dict[Future, str] = {}
        failed: List[SaveFileModel] = []
        for file in file_list:
            try:
                future = self._executor.submit(self._do_upload, file)
            except RuntimeError:
                failed.append(file)
            else:
                submitted[future] = file.name
        if submitted:
            with self._lock:
                self._futures.update(submitted)
        # 线程池已关闭或解释器正在关闭时，逐个手动上传
        self._upload_sync(failed)

    def join(self) -> None:
        while True:
            with self._lock:
                items = list(self._futures.items())
                self._futures.clear()
            if not items:
                return
            wait([f for f, _ in items])
            for future, name in items:
                try:
                    future.result()
                except Exception as e:
                    swanlog.warning(f"Save failed for {name}: {e}")
                if self._progress is not None:
                    self._progress.done()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self.join()
        self._executor.shutdown(wait=True)
        if self._progress is not None:
            self._progress.stop()

    def _upload_sync(self, files: Iterable[SaveFileModel]) -> None:
        for file in files:
            try:
                self._do_upload(file)
            except Exception as e:
                swanlog.warning(f"Save failed for {file.name}: {e}")
            if self._progress is not None:
                self._progress.done()

    def _do_upload(self, file: SaveFileModel) -> None:
        if self._mode == "disabled":
            return
        if not os.path.exists(file.source_path):
            swanlog.warning(f"File not found: {file.source_path}")
            return
        if os.path.getsize(file.source_path) > MAX_FILE_SIZE:
            swanlog.warning(
                f"File '{file.name}' ({file.source_path}) exceeds the size limit "
                f"({MAX_FILE_SIZE // (1024 ** 3)} GB) and will not be uploaded."
            )
            return
        if self._mode in ("local", "offline"):
            self._copy_local(file)
            return
        if self._mode == "cloud":
            self._upload_cloud(file)

    def _copy_local(self, file: SaveFileModel) -> None:
        copy_file(file.source_path, file.target_path)

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
            swanlog.warning(f"Failed to upload {file.name}: {e}")
            self._mark_failed(file, exp_id)

    def _upload_multipart(self, file: SaveFileModel, size: int, exp_id: str) -> None:
        assert self._client is not None
        try:
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
            if upload_id is None:
                raise ValueError("Multipart prepare response is missing uploadId.")
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
        except Exception as e:
            swanlog.warning(f"Failed to upload {file.name}: {e}")
            self._mark_failed(file, exp_id)

    def _mark_failed(self, file: SaveFileModel, exp_id: str) -> None:
        assert self._client is not None
        try:
            complete_upload(
                self._client,
                exp_id,
                [file.complete_request(state=SaveFileState.FAILED)],
            )
        except Exception as e:
            swanlog.warning(f"Failed to report FAILED state for {file.name}: {e}")


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
        self._target_modes: Dict[str, str] = {}
        self._target_modes_lock = threading.Lock()
        self._lock = threading.Lock()
        self._timer = Timer(task=self._poll_task, interval=self._interval)

    def watch(self, files: Union[SaveFileModel, Iterable[SaveFileModel]]) -> None:
        with self._lock:
            for file in _iter_files(files):
                if self._mode == "cloud":
                    self._refresh_live_target(file.name, file)
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
            if self._mode == "cloud":
                self._refresh_live_target(name, current)
            try:
                self._on_change(current)
            except Exception as e:
                swanlog.warning(f"Live save change failed for {name}: {e}")

    def _refresh_live_target(self, name: str, file: SaveFileModel) -> None:
        """刷新 live 模式下的本地镜像，异常时只记录告警，不影响上传流程。"""
        try:
            self._sync_live_target(file)
        except Exception as e:
            swanlog.warning(f"Live save mirror failed for {name}: {e}")

    def _sync_live_target(self, file: SaveFileModel) -> None:
        """根据缓存的同步模式维护 target_path，与 source_path 保持一致。"""
        if same_path(file.source_path, file.target_path):
            return

        mode = self._get_target_mode(file.name)
        if mode == "symlink" and not self._is_symlink_target(file):
            mode = None
            self._clear_target_mode(file.name)
        if mode is None:
            mode = self._resolve_sync_mode(file)

        if mode == "symlink":
            return
        # mode == "copy"
        copy_file(file.source_path, file.target_path)

    def _resolve_sync_mode(self, file: SaveFileModel) -> str:
        """首次确定使用 symlink 还是 copy，并在必要时创建对应的本地镜像。"""
        os.makedirs(os.path.dirname(file.target_path), exist_ok=True)

        if self._is_symlink_target(file):
            self._set_target_mode(file.name, "symlink")
            return "symlink"

        pathlib.Path(file.target_path).unlink(missing_ok=True)
        try:
            rel_target = os.path.relpath(file.source_path, os.path.dirname(file.target_path))
            os.symlink(rel_target, file.target_path)
            self._set_target_mode(file.name, "symlink")
            return "symlink"
        except (NotImplementedError, OSError, ValueError) as e:
            swanlog.debug(
                f"Symlink unavailable for live save target {file.target_path}: {e}. "
                "Falling back to file copy."
            )
            copy_file(file.source_path, file.target_path)
            self._set_target_mode(file.name, "copy")
            return "copy"

    def _get_target_mode(self, name: str) -> Optional[str]:
        with self._target_modes_lock:
            return self._target_modes.get(name)

    def _set_target_mode(self, name: str, mode: str) -> None:
        with self._target_modes_lock:
            self._target_modes[name] = mode

    def _clear_target_mode(self, name: str) -> None:
        with self._target_modes_lock:
            self._target_modes.pop(name, None)

    def _is_symlink_target(self, file: SaveFileModel) -> bool:
        """检查 target_path 是否是一个有效且指向 source_path 的软链接。"""
        if os.path.islink(file.target_path):
            target = os.readlink(file.target_path)
            if not os.path.isabs(target):
                target = os.path.abspath(
                    os.path.join(os.path.dirname(file.target_path), target)
                )
            if same_path(target, file.source_path):
                return True
        return False


__all__ = ["FileUploadManager", "DirWatcher"]
