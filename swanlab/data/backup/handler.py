"""
@author: cunyue
@file: handler.py
@time: 2025/6/4 12:43
@description: 备份处理器，负责日志编解码和写入操作
"""

import os
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional

from swankit.callback import RuntimeInfo, MetricInfo

from swanlab.data.backup.datastore import DataStore
from swanlab.proto.v0 import Experiment, Log, Project, Header, Footer, BaseModel
from swanlab.toolkit import create_time


class BackupHandler:
    """
    备份处理器，负责处理日志备份相关的操作
    """

    BACKUP_FILE = "backup.swanlab"

    def __init__(self):
        super().__init__()
        # 线程执行器
        self.executor: Optional[ThreadPoolExecutor] = None
        # 日志文件写入句柄
        self.f = DataStore()
        # 运行时文件备份目录
        self.files_dir: Optional[str] = None

        # 动态设置包括项目名在内的一些属性，因为在 on_run 之前句柄还未创建，所以需要先缓存，等执行对应的函数的时候再使用
        self.cache_proj_name = None
        self.cache_workspace = None
        self.cache_public = None

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
        self.f.open_for_write(os.path.join(run_dir, self.BACKUP_FILE))
        self.f.write(
            Header.model_validate(
                {
                    "create_time": create_time(),
                    "backup_type": "DEFAULT",
                }
            ).to_record()
        )
        self.backup_proj()
        self.backup_exp(exp_name, description, tags)

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
            self.f.write(log.to_record())
        # 写入结束标志
        footer = Footer.model_validate({"create_time": create_time(), "success": error is None})
        self.f.write(footer.to_record())
        # 关闭日志文件句柄
        self.f.ensure_flushed()
        self.f.close()
        self.f = None

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
        self.f.write(project.to_record())

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
        self.f.write(experiment.to_record())

    def backup(self, data: BaseModel):
        self.f.write(data.to_record())

    def write_runtime_info(self, runtime_info: RuntimeInfo):
        """
        写入运行时信息
        """
        if runtime_info.requirements is not None:
            runtime_info.requirements.write(self.files_dir)
        if runtime_info.metadata is not None:
            runtime_info.metadata.write(self.files_dir)
        if runtime_info.config is not None:
            runtime_info.config.write(self.files_dir)
        if runtime_info.conda is not None:
            runtime_info.conda.write(self.files_dir)

    @staticmethod
    def write_media_buffer(metric_info: MetricInfo):
        """
        写入媒体信息
        """
        if metric_info.metric_buffers is None:
            return
        for i, r in enumerate(metric_info.metric_buffers):
            if r is None:
                continue
            # 组合路径
            path = os.path.join(metric_info.swanlab_media_dir, metric_info.column_info.kid)
            os.makedirs(path, exist_ok=True)
            # 写入数据
            with open(os.path.join(path, metric_info.metric["data"][i]), "wb") as f:
                f.write(r.getvalue())
