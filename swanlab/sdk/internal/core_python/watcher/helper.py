"""
基于 watchdog 的文件监听器辅助定义。

包含 FileEntry、签名计算、回调类型别名、watchdog 事件处理器等支撑组件。
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Optional

from watchdog.events import FileSystemEventHandler

from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.sdk.internal.pkg import fs, safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.watcher import FileWatcher


@dataclass
class FileEntry:
    """已注册文件的追踪信息"""

    name: str  # 相对路径（文件标识）
    source_path: str  # 源文件绝对路径
    target_path: str  # 本地镜像路径（软链接位置）
    policy: Optional[int] = None  # SavePolicy enum int value
    signature: Optional[str] = None  # 上次已知的 mtime+size 签名


@safe.decorator(OSError, level="debug", message=None)
def compute_signature(path: str) -> Optional[str]:
    """基于 mtime + size 计算快速签名，避免每次读文件算 MD5。"""
    st = os.stat(path)
    return f"{st.st_mtime_ns}:{st.st_size}"


OnChangeCallback = Callable[[SaveRecord], None]
"""文件变化回调类型，接收 (SaveRecord)"""


class _Handler(FileSystemEventHandler):
    """watchdog 事件处理器，转发到 FileWatcher 的 debounce 逻辑。"""

    def __init__(self, watcher: "FileWatcher"):
        self._watcher = watcher

    def on_modified(self, event):
        if not event.is_directory:
            self._watcher._schedule_debounce(str(event.src_path))

    def on_created(self, event):
        if not event.is_directory:
            self._watcher._schedule_debounce(str(event.src_path))


def create_save_links(saves: List[SaveRecord], files_dir: Path) -> int:
    """为 SaveRecord 创建软链接并填充 target_path，返回新建链接数量。"""
    count = 0
    for s in saves:
        target = files_dir / s.name
        if not target.exists():
            fs.safe_link(s.source_path, target)
            count += 1
        s.target_path = str(target)
    return count
