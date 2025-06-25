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

from pydantic import BaseModel


class RunStore(BaseModel):

    # ---------------------------------- 项目 ----------------------------------
    # 项目名称
    project: str | None = None
    # 项目所在空间
    workspace: str | None = None
    # 项目可见性
    visibility: bool | None = None
    # ---------------------------------- 实验 ----------------------------------
    # 实验名称
    run_name: str | None = None
    # 实验颜色
    run_colors: list[str] | None = None
    # 实验标签
    tags: list[str] | None = None
    # 实验描述
    description: str | None = None
    # 实验运行 ID
    run_id: str | None = None

    # ---------------------------------- 目录 ----------------------------------
    # 是否为临时目录，标识一些运行时环境
    tmp_dir: bool | None = None
    # 日志存放目录
    swanlog_dir: str | None = None
    # 运行目录
    run_dir: str | None = None

    @property
    def backup_file(self):
        assert os.path.exists(self.run_dir), "Run directory does not exist when accessing backup file."
        return os.path.join(self.run_dir, "backup.swanlab")

    @property
    def log_dir(self):
        assert os.path.exists(self.run_dir), "Run directory does not exist when accessing log directory."
        return os.path.join(self.run_dir, "logs")

    @property
    def console_dir(self):
        assert os.path.exists(self.run_dir), "Run directory does not exist when accessing console directory."
        return os.path.join(self.run_dir, "console")

    @property
    def file_dir(self):
        assert os.path.exists(self.run_dir), "Run directory does not exist when accessing file directory."
        return os.path.join(self.run_dir, "files")

    @property
    def media_dir(self):
        assert os.path.exists(self.run_dir), "Run directory does not exist when accessing media directory."
        return os.path.join(self.run_dir, "media")


run_store = RunStore()


def inside(func):
    """
    检查当前代码是否在 swanlab/data 模块下运行
    如果 swanlab 正在运行测试，则允许在测试中调用此函数
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        frame = inspect.currentframe()
        try:
            caller_module = frame.f_back.f_globals.get('__name__', '')
            if not caller_module.startswith('swanlab.data'):
                if 'PYTEST_VERSION' not in os.environ:
                    raise RuntimeError("This function can only be called from swanlab.data module.")
        finally:
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


__all__ = ["get_run_store", "reset_run_store", "RunStore"]
