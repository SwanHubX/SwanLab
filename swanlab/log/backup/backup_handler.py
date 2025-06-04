"""
@author: cunyue
@file: backup_handler.py
@time: 2025/6/4 12:43
@description: 备份处理器，负责出入日志备份写入操作
"""

import os.path
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional, TextIO

import wrapt
from swankit.callback import ColumnInfo, MetricInfo, RuntimeInfo
from swankit.env import create_time

from swanlab.log.backup.models import Experiment, Log, Project, Column, Runtime, Metric, Header
from swanlab.log.backup.writer import write_media_buffer, write_runtime_info
from swanlab.log.type import LogData
from swanlab.package import get_package_version


def enable_check():
    """
    饰器工厂，根据实例属性决定是否执行函数
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        if getattr(instance, "enable", False):
            return wrapped(*args, **kwargs)
        return None

    return wrapper


def async_io(sync: bool = False):
    """
    类实例异步 IO 方法装饰器，判断是否需要备份
    BackupHandler 实例携带一个线程池，使用此装饰器可以将被装饰的方法放入线程池中执行
    这样能够避免在主线程中执行 IO 密集型操作，提升性能
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        if not getattr(instance, "enable", False):
            return
        if getattr(instance, "f", None) is None:
            return
        # 与 https://github.com/SwanHubX/SwanLab/issues/889 相同的问题
        # 在回调中线程池已经关闭，我们需要在主线程中执行
        if sync:
            return wrapped(*args, **kwargs)

        executor: Optional[ThreadPoolExecutor] = getattr(instance, "executor")
        if executor is None or executor._shutdown:
            return wrapped(*args, **kwargs)
        executor.submit(wrapped, *args, **kwargs)

    return wrapper


class BackupHandler:
    """
    备份处理器，负责处理日志备份相关的操作
    """

    def __init__(self, enable: bool = True, backup_type: str = "DEFAULT", save_media: bool = True):
        super().__init__()
        # 是否启用备份
        self.enable = enable
        self.backup_type = backup_type
        # 线程执行器
        self.executor: Optional[ThreadPoolExecutor] = None
        # 日志文件写入句柄
        self.f: Optional[TextIO] = None
        # 运行时文件备份目录
        self.files_dir: Optional[str] = None
        self.save_media: bool = save_media

        # 动态设置包括项目名在内的一些属性，因为在 on_run 之前句柄还未创建，所以需要先缓存，等执行对应的函数的时候再使用
        self.cache_proj_name = None
        self.cache_workspace = None
        self.cache_public = None

    @enable_check()
    def start(self, run_dir: str, files_dir: str, exp_name: str, description: str, tags: List[str]):
        """
        开启备份处理器，创建日志文件句柄
        此函数的功能包括：
        1. 创建线程执行器
        2. 创建日志文件句柄
        3. 在日志文件头写入当前备份类型和一些元信息
        4. 写入项目、实验备份
        """
        self.files_dir = files_dir
        # 创建线程池执行器，每次只会有一个线程在执行，这样的设计原因为：
        # 1. io 操作不会特别耗时
        # 2. 避免多线程写入同一文件导致数据混乱
        # 3. 部分用户会将 swanlog 文件夹挂载在 NAS 等对写入并发有限制的存储设备上
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.f = open(os.path.join(run_dir, "backup.swanlab"), "a", encoding="utf-8")
        self.f.write(
            Header.model_validate(
                {
                    "create_time": create_time(),
                    "version": get_package_version(),
                    "backup_type": self.backup_type,
                }
            ).to_backup()
            + "\n"
        )
        self.backup_proj()
        self.backup_exp(exp_name, description, tags)

    @enable_check()
    def stop(self, epoch: int, error: str = None):
        """
        停止备份处理器
        :param epoch: int, 日志行数
        :param error: str, 如果有错误信息，则在日志中记录
        """
        # 同步停止
        self.executor.shutdown(wait=True)
        # 如果有错误信息则在日志中记录
        if error is not None:
            log = Log.model_validate({"level": "ERROR", "message": error, "create_time": create_time(), "epoch": epoch})
            self.f.write(log.to_backup() + "\n")
        # 关闭日志文件句柄
        self.f.close()
        self.f = None

    @async_io()
    def backup_terminal(self, log_data: LogData):
        """
        备份终端输出
        """
        logs = Log.from_log_data(log_data)
        for log in logs:
            self.f.write(log.to_backup() + "\n")

    @async_io()
    def backup_proj(self):
        """
        备份项目信息
        """
        project = Project.model_validate(
            {
                "name": self.cache_proj_name,
                "workspace": self.cache_workspace,
                "public": self.cache_public,
            }
        )
        self.f.write(project.to_backup() + "\n")

    @async_io()
    def backup_exp(self, exp_name: str, description: str, tags: List[str]):
        """
        备份实验信息
        """
        experiment = Experiment.model_validate(
            {
                "name": exp_name,
                "description": description,
                "tags": tags,
            }
        )
        self.f.write(experiment.to_backup() + "\n")

    @async_io()
    def backup_column(self, column_info: ColumnInfo):
        """
        备份指标列信息
        """
        column = Column.from_column_info(column_info)
        self.f.write(column.to_backup() + "\n")

    @async_io()
    def backup_runtime(self, runtime_info: RuntimeInfo):
        """
        备份运行时信息
        """
        runtime = Runtime.from_runtime_info(runtime_info)
        self.f.write(runtime.to_backup() + "\n")
        write_runtime_info(self.files_dir, runtime_info)

    @async_io()
    def backup_metric(self, metric_info: MetricInfo):
        """
        备份指标信息
        """
        metric = Metric.from_metric_info(metric_info)
        self.f.write(metric.to_backup() + "\n")
        if self.save_media:
            write_media_buffer(metric_info)


def backup(method: str):
    """
    备份装饰器，用于在方法执行前进行备份操作
    """

    @wrapt.decorator
    def wrapper(wrapped, obj, args, kwargs):
        # 执行备份操作
        backup_obj = getattr(obj, "backup")
        getattr(backup_obj, f"backup_{method}")(*args, **kwargs)
        # 执行原方法
        wrapped(*args, **kwargs)

    return wrapper
