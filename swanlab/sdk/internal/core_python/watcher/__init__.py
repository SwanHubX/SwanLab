"""
基于 watchdog 的文件监听器，采用 trailing debounce 策略。

监听 swanlog/{run_id}/files/ 目录下的文件变化，
文件稳定（停止写入 debounce_delay 秒）后触发 on_change 回调。
"""

import os
import threading
from pathlib import Path
from typing import Dict, List, Optional

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
        self._registered: Dict[str, FileEntry] = {}  # abs_path → entry
        self._lock = threading.Lock()
        self._started = False

    def watch(self, dir_path: str, file_paths: List[str], policies: Optional[List[int]] = None) -> None:
        """注册并开始监听指定目录下的文件。

        :param dir_path: 监听目录的绝对路径
        :param file_paths: 相对于 dir_path 的文件路径列表
        :param policies: SavePolicy enum values corresponding to file_paths
        """
        dir_abs = str(Path(dir_path).resolve())

        # 注册文件，计算初始签名
        with self._lock:
            for idx, rel in enumerate(file_paths):
                abs_path = str(Path(dir_abs) / rel)
                if abs_path in self._registered:
                    continue
                source_path = self._resolve_source(abs_path)
                entry = FileEntry(
                    name=rel,
                    source_path=source_path,
                    target_path=abs_path,
                    policy=policies[idx] if policies and idx < len(policies) else None,
                    signature=compute_signature(abs_path),
                )
                self._registered[abs_path] = entry

        # 启动 watchdog（只需一次）
        if not self._started:
            self._observer.schedule(_Handler(self), dir_abs, recursive=False)
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
        """定时器到期后执行：计算签名 → 对比 → 触发回调。"""
        with self._lock:
            self._timers.pop(path, None)
            entry = self._registered.get(path)
            if entry is None:
                return

        # 文件被删除则移除注册
        new_sig = compute_signature(path)
        if new_sig is None:
            with self._lock:
                self._registered.pop(path, None)
            return

        # 签名未变则忽略
        if new_sig == entry.signature:
            return

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
        """对 policy=SAVE_POLICY_LIVE 的记录注册文件监听。"""
        live_files = [s for s in save_records if s.policy == SavePolicy.SAVE_POLICY_LIVE]
        if not live_files:
            return
        self.watch(str(files_dir), [s.name for s in live_files], [s.policy for s in live_files])

    def stop(self) -> None:
        """停止监听，释放资源。"""
        with self._lock:
            for timer in self._timers.values():
                timer.cancel()
            self._timers.clear()

        if self._started:
            self._observer.stop()
            with safe.block(message="FileWatcher observer join error"):
                self._observer.join(timeout=5)
            self._started = False


__all__ = ["create_save_links", "FileWatcher"]
