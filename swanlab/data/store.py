"""
@author: cunyue
@file: store.py
@time: 2025/6/20 17:00
@description: 运行时配置
NOTE: 只允许在 swanlab/data 模块下访问，其他地方不允许访问
"""

import functools
import inspect
import os.path
import time
from typing import Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel

# 定义远程指标类型
RemoteMetric = Dict[str, Tuple[str, str, Optional[dict], Optional[int]]]

_RUN_DIR_VERIFY_RETRIES = 20
_RUN_DIR_VERIFY_INTERVAL = 0.05


class RunStore(BaseModel):

    # ---------------------------------- 项目 ----------------------------------
    # 项目名称
    project: Optional[str] = None
    # 项目所在空间
    workspace: Optional[str] = None
    # 项目可见性
    visibility: Optional[bool] = None
    # ---------------------------------- 实验 ----------------------------------
    # 实验模式
    resume: Literal['must', 'allow', 'never'] = 'never'
    # 实验名称
    run_name: Optional[str] = None
    # 实验颜色
    run_colors: Optional[Tuple[str, str]] = None
    # 实验标签
    tags: Optional[List[str]] = None
    # 任务类型
    job_type: Optional[str] = None
    # 实验组
    group: Optional[str] = None
    # 实验描述
    description: Optional[str] = None
    # 实验运行 ID
    run_id: Optional[str] = None
    # 当前实验是否为新实验
    new: Optional[bool] = None
    # 恢复实验时，云端实验的 config 设置
    config: Optional[dict] = None
    # 恢复实验时，云端实验的指标数据，key -> (column_type, column_class, error, latest step)
    metrics: Optional[RemoteMetric] = None
    # 恢复实验时，云端实验的日志条数
    log_epoch: Optional[int] = None
    # ---------------------------------- 目录 ----------------------------------
    # 是否为临时目录，标识一些运行时环境
    tmp_dir: Optional[bool] = None
    # 日志存放目录
    swanlog_dir: Optional[str] = None
    # 运行目录
    run_dir: Optional[str] = None

    def _ensure_run_dir(self) -> str:
        """
        统一的 run_dir 校验入口。
        1. run_dir 不能为 None（未初始化）
        2. run_dir 必须存在（被删除或路径错误），短暂不可见时等待重试
        首次校验通过后缓存结果，避免频繁 I/O。
        """
        run_dir = self.run_dir
        if run_dir is None:
            raise RuntimeError(
                "Run directory has not been initialized. "
                "This means swanlab.init() has not been called, or the run has already been finished."
            )
        if not self._run_dir_verified:
            self._wait_for_run_dir(run_dir)
            self._run_dir_verified = True
        return run_dir

    def _wait_for_run_dir(self, run_dir: str) -> None:
        for i in range(_RUN_DIR_VERIFY_RETRIES + 1):
            if os.path.isdir(run_dir):
                return
            if i < _RUN_DIR_VERIFY_RETRIES:
                time.sleep(_RUN_DIR_VERIFY_INTERVAL)

        raise FileNotFoundError(
            f"Run directory does not exist: {run_dir}. "
            f"The experiment data may have been moved, deleted, or the path is incorrect."
        )

    # Pydantic private attribute, not included in schema or serialization.
    _run_dir_verified: bool = False

    @property
    def backup_file(self) -> str:
        return os.path.join(self._ensure_run_dir(), "backup.swanlab")

    @property
    def log_dir(self) -> str:
        return os.path.join(self._ensure_run_dir(), "logs")

    @property
    def console_dir(self) -> str:
        return os.path.join(self._ensure_run_dir(), "console")

    @property
    def file_dir(self) -> str:
        return os.path.join(self._ensure_run_dir(), "files")

    @property
    def media_dir(self) -> str:
        return os.path.join(self._ensure_run_dir(), "media")


run_store = RunStore()


def inside(func):
    """
    检查当前代码是否在 swanlab/data 模块下运行
    如果 swanlab 正在运行测试，则允许在测试中调用此函数
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        frame = inspect.currentframe()
        caller_frame = frame.f_back if frame is not None else None
        try:
            caller_module = caller_frame.f_globals.get('__name__', '') if caller_frame is not None else ''
            if not caller_module.startswith('swanlab.data') and 'PYTEST_VERSION' not in os.environ:
                raise RuntimeError("This function can only be called from swanlab.data module.")
        finally:
            del caller_frame
            del frame
        return func(*args, **kwargs)

    return wrapper


@inside
def get_run_store():
    """
    此模块只允许在 swanlab/data 模块下访问
    为了提高性能，建议尽量减少对此函数的调用次数
    """
    global run_store
    return run_store


@inside
def reset_run_store():
    global run_store
    run_store = RunStore()
    return None


__all__ = ["get_run_store", "reset_run_store", "RunStore", "RemoteMetric"]
