"""
@author: cunyue
@file: sync_utils.py
@time: 2025/7/21 17:20
@description: swanlab sync 函数的工具函数
"""

import threading
from typing import Optional, Union, Literal

from rich.progress import (
    Progress,
    BarColumn,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
    Task,
)

from ..data.store import RunStore
from ..proto.v0 import Project, Experiment


class SyncProgress:
    """
    同步进度条上下文管理器，处理进度条展示
    """

    def __init__(self):
        self._pbar = Progress(
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
            MofNCompleteColumn(),
        )
        self._task = None
        self._lock = threading.Lock()

    def __enter__(self):
        self._pbar.start()
        return self

    @property
    def task(self) -> Task:
        t = getattr(self._pbar, '_tasks', {}).get(self._task)
        if t is None:
            raise RuntimeError("No task found.")
        return t

    def __exit__(self, exc_type, exc_val, exc_tb):
        # 如果没有错误，确保进度条到达100%
        if exc_type is None and self._task is not None:
            with self._lock:
                if self.task.completed < self.task.total:
                    self._pbar.update(self._task, completed=self.task.total)
        self._pbar.stop()

    def set_total(self, total: int):
        """
        设置进度条总数并开始任务
        """
        self._task = self._pbar.add_task("Syncing data...", total=total)

    def update(self, uploaded: int):
        """
        进度回调函数
        :param uploaded: 本次已上传的数量
        """
        print("uploaded:", uploaded)
        if self._task is not None and uploaded > 0:
            with self._lock:
                self._pbar.update(self._task, advance=uploaded)


def set_run_store(
    run_store: RunStore,
    proj: Project,
    exp: Experiment,
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    id: Union[str, Literal['auto', 'new']] = "new",
):
    """
    在同步之前设置运行存储的相关信息。
    :param run_store: RunStore 实例
    :param proj: 项目信息
    :param exp: 实验信息
    :param project: 指定的项目名称，如果为 None 则使用 proj.name
    :param workspace: 指定的工作空间名称，如果为 None 则使用 proj.workspace
    :param id: 实验 ID，可以是 'new', 'auto' 或具体的 ID 字符串 本函数不做字符串格式检查
    """
    # 设置项目信息
    run_store.project = project or proj.name
    run_store.workspace = workspace or proj.workspace
    run_store.visibility = proj.public
    # 设置实验信息
    run_store.run_name = exp.name
    run_store.description = exp.description
    run_store.tags = exp.tags
    run_store.run_colors = (exp.colors[0], exp.colors[1])
    # 设置实验 id 和 resume 模式
    # a. id 为 new，则 resume 为 never, id 为 None
    # b. id 为 auto，则 resume 为 allow, id 为 exp.id
    # c. 其他情况，则 resume 为 must, id 为 id
    if id == "new":
        run_store.resume = 'never'
        run_store.run_id = None
    elif id == "auto":
        run_store.resume = 'allow'
        run_store.run_id = exp.id
    else:
        run_store.resume = 'must'
        run_store.run_id = id
