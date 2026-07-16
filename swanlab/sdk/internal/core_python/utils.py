"""
@author: cunyue
@file: utils.py
@time: 2026/5/19 15:55
@description: Core 服务共享工具函数
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional

from swanlab.proto.swanlab.run.v1.run_pb2 import StartRecord
from swanlab.sdk.internal.core_python.api.experiment import (
    create_or_resume_experiment as api_create_or_resume_experiment,
)
from swanlab.sdk.internal.core_python.api.project import get_or_create_project
from swanlab.sdk.internal.pkg import adapter
from swanlab.sdk.typings.core_python.api.experiment import InitExperimentType
from swanlab.sdk.typings.core_python.api.project import ProjectType
from swanlab.utils.experiment import generate_color, generate_name


@dataclass(frozen=True)
class PrepareExperimentStartResult:
    username: str
    project: str
    project_data: ProjectType
    experiment: InitExperimentType
    new_experiment: bool
    name: str
    color: str


def prepare_experiment_start(record: StartRecord) -> PrepareExperimentStartResult:
    """
    创建或恢复运行对应的项目和实验
    """
    project_data = get_or_create_project(
        username=record.workspace,
        name=record.project,
        public=record.public,
    )
    username, project = project_data["group"]["username"], project_data["name"]
    history_experiment_count = project_data["_count"]["experiments"]
    name = record.name or generate_name(history_experiment_count)
    color = record.color or generate_color(history_experiment_count)
    resume = adapter.resume[record.resume]

    experiment, new_experiment = api_create_or_resume_experiment(
        username,
        project,
        name=name,
        resume=resume,
        run_id=record.id,
        color=color,
        description=record.description,
        job_type=record.job_type,
        group=record.group,
        tags=list(record.tags),
        created_at=record.started_at,
    )
    assert experiment.get("name"), "create_or_resume_experiment() returned an experiment without a name."

    return PrepareExperimentStartResult(
        username=username,
        project=project,
        project_data=project_data,
        experiment=experiment,
        new_experiment=new_experiment,
        name=name,
        color=color,
    )


def generate_run_online_path(result: PrepareExperimentStartResult):
    # 获取path，/:username/:project_name/:run_id
    run_id = result.experiment.get("slug", "") or result.experiment["cuid"]
    # /:username/:project_name
    project_path = result.project_data["path"]
    return f"{project_path}/{run_id}"


def get_buffer_size(buffer: Any) -> int:
    """获取文件类对象的字节大小，支持 str/Path（文件路径）、BytesIO、文件句柄等。

    对于有 fileno 的真实文件优先走 fstat（无副作用）；
    对于内存 buffer 走 getbuffer / getvalue / __len__；
    最后用 seek+tell 回退方案。无法确定大小时抛 TypeError。
    """
    # 文件路径：直接 stat，不打开文件
    if isinstance(buffer, (str, Path)):
        return os.path.getsize(buffer)
    # BytesIO / StringIO
    if hasattr(buffer, "getbuffer"):
        return buffer.getbuffer().nbytes
    if hasattr(buffer, "getvalue"):
        return len(buffer.getvalue())
    # 通用有 __len__ 的对象
    if hasattr(buffer, "__len__"):
        return len(buffer)
    # 某些库（如 requests）用 .len 属性
    if hasattr(buffer, "len"):
        return getattr(buffer, "len")
    # 真实文件句柄：fstat 不会改变文件指针
    try:
        return os.fstat(buffer.fileno()).st_size
    except (AttributeError, ValueError, OSError):
        pass
    # 最后手段：seek 到末尾再恢复位置（有副作用，但别无选择）
    if hasattr(buffer, "seek") and hasattr(buffer, "tell"):
        curr = buffer.tell()
        buffer.seek(0, os.SEEK_END)
        size = buffer.tell()
        buffer.seek(curr)
        return size
    raise TypeError("Object has no len")


class ProgressFileWrapper:
    """包装文件对象，在每次 read() 后回调汇报已读字节数。

    用于 requests session.put(data=wrapper) 时自动追踪上传字节进度。
    on_read 回调接收当前累计已读字节数（clamp 到 total_size 以内）。

    支持 offset + size 切片：传入 offset 后仅暴露文件中 [offset, offset+size) 的区间，
    read/tell/seek 均相对于该切片，不再需要额外的 FileReader 包装层。
    """

    def __init__(
        self,
        file_obj: Any,
        on_read: Optional[Callable[[int], None]] = None,
        total_size: Optional[int] = None,
        offset: Optional[int] = None,
        size: Optional[int] = None,
    ):
        self._file_obj = file_obj
        self._on_read = on_read
        self._slice_size = max(size, 0) if size is not None else None
        self._base_offset = max(offset, 0) if offset is not None else None
        self._total_size = self._slice_size if total_size is None else total_size
        if self._base_offset is not None:
            self._file_obj.seek(self._base_offset)
            self._current = 0
        else:
            self._current = self._clamp_position(self.tell())

    def read(self, size: int = -1) -> bytes:
        # 有切片时限制读取范围
        if self._slice_size is not None:
            remaining = self._slice_size - self._current
            if remaining <= 0:
                return b""
            if size is None or size < 0 or size > remaining:
                size = remaining
        chunk = self._file_obj.read(size)
        if chunk:
            self._current = self._clamp_position(self._current + len(chunk))
            if self._on_read:
                self._on_read(self._current)
        return chunk

    def seek(self, offset: int, whence: int = 0) -> int:
        if self._slice_size is not None:
            # 切片模式下 seek 相对于切片内部
            if whence == 0:
                position = offset
            elif whence == 1:
                position = self._current + offset
            elif whence == 2:
                position = self._slice_size + offset
            else:
                raise ValueError(f"Invalid whence: {whence}")
            self._current = min(max(int(position), 0), self._slice_size)
            underlying = (self._base_offset or 0) + self._current
            self._file_obj.seek(underlying)
            return self._current
        position = self._file_obj.seek(offset, whence)
        self._current = self._clamp_position(position)
        return position

    def tell(self) -> int:
        if self._slice_size is not None:
            return self._current
        return self._file_obj.tell()

    def _clamp_position(self, position: int) -> int:
        position = max(0, int(position))
        if self._total_size is not None:
            return min(position, self._total_size)
        return position

    def __getattr__(self, name: str) -> Any:
        return getattr(self._file_obj, name)

    def __len__(self) -> int:
        if self._slice_size is not None:
            return self._slice_size
        return get_buffer_size(self._file_obj)
