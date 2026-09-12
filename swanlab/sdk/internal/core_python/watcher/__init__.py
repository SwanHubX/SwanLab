"""
基于 watchdog 的文件监听器，采用 trailing debounce 策略。

两种注册模式：
  1. 镜像模式（默认）：监听 swanlog/{run_id}/files/ 下的软链接镜像；
  2. direct-source（skip_store 下）：没有本地镜像目录，直接监听用户源文件所在目录，
     事件路径与注册的源文件绝对路径精确匹配，同目录其他文件的变化被忽略。

文件稳定（停止写入 debounce_delay 秒）后触发 on_change 回调。
"""

import os
import threading
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from watchdog.observers import Observer

from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.sdk.internal.pkg import safe

from .helper import FileEntry, OnChangeCallback, _Handler, compute_signature, create_save_links


class FileWatcher:
    """基于 watchdog + trailing debounce 的文件监听器。

    参数:
        on_change: 文件变化后的回调，接收 (abs_path, SaveRecord)
        debounce_delay: 文件停止变化后等待多少秒再触发回调，默认 1.0
    """

    def __init__(self, on_change: OnChangeCallback, debounce_delay: float = 1.0):
        self._on_change = on_change
        self._debounce_delay = debounce_delay
        self._observer = Observer()
        self._timers: Dict[str, threading.Timer] = {}
        self._registered: Dict[str, List[FileEntry]] = {}  # 事件路径 → entries（同源多 name 时为多条）
        self._scheduled_dirs: Set[str] = set()  # 已调度监听的目录绝对路径
        self._lock = threading.Lock()
        self._started = False

    def watch(self, dir_path: str, file_paths: List[str], policies: Optional[List[int]] = None) -> None:
        """镜像模式：注册并开始监听指定目录下的文件。

        :param dir_path: 监听目录的绝对路径
        :param file_paths: 相对于 dir_path 的文件路径列表
        :param policies: SavePolicy enum values corresponding to file_paths
        """
        dir_abs = str(Path(dir_path).resolve())

        # 注册文件，计算初始签名
        with self._lock:
            for idx, rel in enumerate(file_paths):
                abs_path = str(Path(dir_abs) / rel)
                source_path = self._resolve_source(abs_path)
                entry = FileEntry(
                    name=rel,
                    source_path=source_path,
                    target_path=abs_path,
                    policy=policies[idx] if policies and idx < len(policies) else None,
                    signature=compute_signature(abs_path),
                )
                self._register(abs_path, entry)

        self._ensure_scheduled(dir_abs)

    def watch_sources(self, saves: List[SaveRecord]) -> None:
        """direct-source 模式：不依赖本地镜像目录，直接监听 source_path 所在目录。

        - 按 source_path.parent 分组，一个目录只 schedule 一次；
        - _registered 以源文件绝对路径为 key，事件精确匹配（同目录其他文件被忽略）；
        - 同一源文件保存为多个 name 时一对多注册，变化时为每个 name 各触发一次回调。
        """
        # 1. 按源文件父目录分组
        groups: Dict[str, List[Tuple[str, SaveRecord]]] = {}
        for save in saves:
            if not save.source_path:
                continue
            source_abs = str(Path(save.source_path).resolve())
            groups.setdefault(str(Path(source_abs).parent), []).append((source_abs, save))
        # 2. 逐目录注册并调度监听
        for dir_abs, group in groups.items():
            with self._lock:
                for source_abs, save in group:
                    self._register(
                        source_abs,
                        FileEntry(
                            name=save.name,
                            source_path=source_abs,
                            target_path="",  # direct-source：无本地镜像路径
                            policy=save.policy,
                            signature=compute_signature(source_abs),
                        ),
                    )
            self._ensure_scheduled(dir_abs)

    def _register(self, event_path: str, entry: FileEntry) -> None:
        """登记一条监听（调用方持锁）。同一事件路径可挂多个 name，按 (name, source_path) 幂等。"""
        entries = self._registered.setdefault(event_path, [])
        if any(e.name == entry.name and e.source_path == entry.source_path for e in entries):
            return
        entries.append(entry)

    def _ensure_scheduled(self, dir_abs: str) -> None:
        """确保目录已被 watchdog 监听；首次调度时启动 observer 线程。"""
        if dir_abs in self._scheduled_dirs:
            return
        self._scheduled_dirs.add(dir_abs)
        self._observer.schedule(_Handler(self), dir_abs, recursive=False)
        if not self._started:
            self._observer.start()
            self._started = True

    def _resolve_source(self, target_path: str) -> str:
        """如果 target_path 是软链接，解析到源文件；否则返回自身。"""
        with safe.block(OSError, level="debug", message=None):
            if os.path.islink(target_path):
                return os.path.realpath(target_path)
        return target_path

    def _schedule_debounce(self, path: str) -> None:
        """为 path 启动/重置 trailing debounce 定时器。"""
        with self._lock:
            if path not in self._registered:
                return
            # 取消已有定时器
            old = self._timers.pop(path, None)
            if old is not None:
                old.cancel()
            # 创建新定时器
            timer = threading.Timer(self._debounce_delay, self._process_change, args=(path,))
            self._timers[path] = timer
        timer.start()

    def _process_change(self, path: str) -> None:
        """定时器到期后执行：计算签名 → 对比 → 触发回调。同一路径的所有条目各回调一次。"""
        with self._lock:
            self._timers.pop(path, None)
            entries = self._registered.get(path)
            if not entries:
                return
            entries = list(entries)

        # 文件被删除则移除注册
        new_sig = compute_signature(path)
        if new_sig is None:
            with self._lock:
                self._registered.pop(path, None)
            return

        for entry in entries:
            # 签名未变则忽略
            if new_sig == entry.signature:
                continue
            # 签名变化，更新并触发回调
            entry.signature = new_sig
            record = SaveRecord(
                name=entry.name,
                source_path=entry.source_path,
                target_path=entry.target_path,
            )
            if entry.policy is not None:
                record.policy = entry.policy  # type: ignore[assignment]
            with safe.block(message=f"FileWatcher on_change callback error for {path}"):
                self._on_change(record)

    def register_live_watches(self, save_records: List[SaveRecord], files_dir: Path) -> None:
        """镜像模式：对 policy=SAVE_POLICY_LIVE 的记录注册文件监听。"""
        live_files = [s for s in save_records if s.policy == SavePolicy.SAVE_POLICY_LIVE]
        if not live_files:
            return
        self.watch(str(files_dir), [s.name for s in live_files], [s.policy for s in live_files])

    def register_source_watches(self, save_records: List[SaveRecord]) -> None:
        """direct-source 模式（skip_store 下）：对 policy=SAVE_POLICY_LIVE 的记录直接监听源文件。"""
        live_files = [s for s in save_records if s.policy == SavePolicy.SAVE_POLICY_LIVE]
        if not live_files:
            return
        self.watch_sources(live_files)

    def stop(self) -> None:
        """停止监听，释放资源。"""
        with self._lock:
            for timer in self._timers.values():
                timer.cancel()
            self._timers.clear()
            self._scheduled_dirs.clear()

        if self._started:
            self._observer.stop()
            with safe.block(message="FileWatcher observer join error"):
                self._observer.join(timeout=5)
            self._started = False


__all__ = ["create_save_links", "FileWatcher"]
