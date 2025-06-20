"""
@author: cunyue
@file: handler.py
@time: 2025/6/4 12:43
@description: 备份处理器，负责日志编解码和写入操作
"""

import os
from concurrent.futures import ThreadPoolExecutor

from swankit.callback import RuntimeInfo, MetricInfo

from swanlab.data.backup.datastore import DataStore
from swanlab.data.store import run_store
from swanlab.proto.v0 import Log, Project, Header, Footer, BaseModel, Experiment
from swanlab.toolkit import create_time


class BackupHandler:
    """
    备份处理器，负责处理日志备份相关的操作
    """

    BACKUP_FILE = "backup.swanlab"

    def __init__(self):
        super().__init__()
        # 创建线程池执行器，每次只会有一个线程在执行，这样的设计原因为：
        # 1. io 操作不会特别耗时
        # 2. 避免多线程写入同一文件导致数据混乱
        # 3. 部分用户会将 swanlog 文件夹挂载在 NAS 等对写入并发有限制的存储设备上
        self.executor = ThreadPoolExecutor(max_workers=1)
        # 日志文件写入句柄
        self.f = DataStore()

    def start(self):
        """
        开启备份处理器，创建日志文件句柄
        此函数的功能包括：
        1. 创建线程执行器
        2. 创建日志文件句柄
        3. 在日志文件头写入当前备份类型和一些元信息
        4. 写入项目、实验备份
        """
        assert os.path.exists(run_store.file_dir)
        self.f.open_for_write(run_store.backup_file)
        self.f.write(
            Header.model_validate(
                {
                    "create_time": create_time(),
                    "backup_type": "DEFAULT",
                }
            ).to_record()
        )
        self.f.write(
            Project.model_validate(
                {
                    "name": run_store.run_name,
                    "workspace": run_store.workspace,
                    "public": run_store.visibility,
                }
            ).to_record()
        )
        self.f.write(
            Experiment.model_validate(
                {
                    "name": run_store.run_name,
                    "description": run_store.description,
                    "tags": run_store.tags,
                }
            ).to_record()
        )

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

    def backup(self, data: BaseModel):
        self.f.write(data.to_record())

    @staticmethod
    def write_runtime_info(runtime_info: RuntimeInfo):
        """
        写入运行时信息
        """
        if runtime_info.requirements is not None:
            runtime_info.requirements.write(run_store.file_dir)
        if runtime_info.metadata is not None:
            runtime_info.metadata.write(run_store.file_dir)
        if runtime_info.config is not None:
            runtime_info.config.write(run_store.file_dir)
        if runtime_info.conda is not None:
            runtime_info.conda.write(run_store.file_dir)

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
