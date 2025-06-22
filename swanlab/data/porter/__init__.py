"""
@author: cunyue
@file: __init__.py.py
@time: 2025/6/22 19:46
@description: 数据搬运工，在 callbacker 中被实例化，支持：
1. 数据上报（可选）：将数据上报到 SwanLab 服务器
2. 数据备份（必须）：将数据备份到本地文件系统
3. 读取本地备份并上报：读取本地备份文件并上报到 SwanLab 服务器

其中数据上报可选择使用 go 后端或者 python 线程
SwanLab 不生产数据，我们只是 AI 训练的搬运工，祝愿所有训练者跑出预期的 loss :)
"""

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Literal

import wrapt
from swankit.callback import MetricInfo, ColumnInfo, RuntimeInfo

from swanlab.core_python import get_client
from swanlab.core_python.uploader.thread import ThreadPool, UploadType
from swanlab.data.store import RunStore, get_run_store
from swanlab.log.type import LogData
from swanlab.proto.v0 import Log, Header, Project, Experiment, Column, Metric, BaseModel, Runtime, Footer
from swanlab.toolkit import create_time
from .datastore import DataStore


def traced():
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        if instance is not None and not getattr(instance, "_traced", False):
            raise RuntimeError("ProtoTransfer has not been started for tracing. Call start_for_trace() first.")
        return wrapped(*args, **kwargs)

    return wrapper


def once():
    called = False

    @wrapt.decorator
    def wrapper(wrapped, _, args, kwargs):
        nonlocal called
        if not called:
            called = True
            return wrapped(*args, **kwargs)
        raise RuntimeError('This method can only be called once.')

    return wrapper


def backup():
    """
    标记需要备份的函数，并将函数的返回结果存储到备份文件中
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        result: BaseModel | list[BaseModel] = wrapped(*args, **kwargs)
        if result is not None:
            f = getattr(instance, "_f")
            if isinstance(result, list):
                [f.write(item.to_record()) for item in result]
            else:
                f.write(result.to_record())
        return result

    return wrapper


def async_io():
    """
    标记需要在其他线程处理的函数
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        executor: ThreadPoolExecutor = getattr(instance, "_executor", None)
        # 如果没有设置 executor，则直接执行函数
        if executor is None:
            return wrapped(*args, **kwargs)
        # 与 https://github.com/SwanHubX/SwanLab/issues/889 相同的问题
        # 在回调中线程池已经关闭，我们需要在主线程中执行
        if executor._shutdown:
            return wrapped(*args, **kwargs)
        # 否则将函数提交到线程池中执行
        executor.submit(wrapped, *args, **kwargs)

    return wrapper


class DataPorter:
    _instance: Optional['DataPorter'] = None
    _run_store: Optional[RunStore] = None

    def __init__(self):
        # 上传线程池，可以选择不开启
        self._pool: Optional[ThreadPool] = None
        # 数据存储句柄
        self._f = DataStore()
        # 是否开启日志跟踪
        self._traced = False
        # 是否开启异步
        self._executor: Optional[ThreadPoolExecutor] = None

    @once()
    def open_for_trace(self, sync: bool = False, backend: Literal['go', 'python', 'none'] = 'python'):
        """
        开启日志传输器以实验日志跟踪，创建数据存储句柄
        :param sync: 是否同步上传数据，默认为 False，表示异步上传
        :param backend: 上传后端类型，默认为 'python'，可选 'go' 或 'none', none 表示不开启上传
        """
        self._f.open_for_write(self._run_store.backup_file)
        if backend == 'python':
            self._pool = ThreadPool()
        elif backend == 'go':
            raise NotImplementedError("swanlab-core is not ready yet.")
        elif backend == 'none':
            pass
        else:
            raise ValueError(f"Unsupported backend: {backend}")

        if not sync:
            self._executor = ThreadPoolExecutor(max_workers=1)

        # 写入备份文件头
        self._f.write(
            Header.model_validate(
                {
                    "create_time": create_time(),
                    "backup_type": "DEFAULT",
                }
            ).to_record()
        )
        run_name = self._run_store.run_name
        workspace = self._run_store.workspace
        visibility = self._run_store.visibility
        description = self._run_store.description
        tags = self._run_store.tags
        self._f.write(
            Project.model_validate(
                {
                    "name": run_name,
                    "workspace": workspace,
                    "public": visibility,
                }
            ).to_record()
        )
        self._f.write(
            Experiment.model_validate(
                {
                    "name": run_name,
                    "description": description,
                    "tags": tags,
                }
            ).to_record()
        )
        self._traced = True

    def _publish(self, *args, **kwargs):
        """
        发布数据到上传线程池，如果未开启线程池，则不执行任何操作
        """
        if self._pool is not None:
            self._pool.queue.put(*args, **kwargs)

    @backup()
    @traced()
    @async_io()
    def trace_column(self, data: ColumnInfo) -> BaseModel:
        """
        追踪列数据
        """
        column = Column.from_column_info(data)
        self._publish((UploadType.COLUMN, [column.to_column_model()]))
        return column

    @backup()
    @traced()
    @async_io()
    def trace_metric(self, data: MetricInfo) -> BaseModel:
        """
        追踪指标数据
        """
        assert data.error is None, "MetricInfo must have an error field"
        if data.column_info.chart_type == data.column_info.chart_type.LINE:
            # 标量
            scalar = Metric.from_metric_info(data)
            self._publish((UploadType.SCALAR_METRIC, [scalar.to_scalar_model()]))
            return scalar
        else:
            # 媒体
            media = Metric.from_metric_info(data)
            # 写入媒体信息
            if data.metric_buffers is not None:
                for i, r in enumerate(data.metric_buffers):
                    if r is None:
                        continue
                    # 组合路径
                    path = os.path.join(data.swanlab_media_dir, data.column_info.kid)
                    os.makedirs(path, exist_ok=True)
                    # 写入数据
                    with open(os.path.join(path, data.metric["data"][i]), "wb") as f:
                        f.write(r.getvalue())
            self._publish((UploadType.MEDIA_METRIC, [media.to_media_model(data.swanlab_media_dir)]))
            return media

    @backup()
    @traced()
    @async_io()
    def trace_runtime_info(self, data: RuntimeInfo) -> BaseModel:
        """
        追踪运行时信息
        """
        runtime = Runtime.from_runtime_info(data)
        # 写入运行时信息到数据存储
        file_dir = self._run_store.file_dir
        if data.requirements is not None:
            data.requirements.write(file_dir)
        if data.metadata is not None:
            data.metadata.write(file_dir)
        if data.config is not None:
            data.config.write(file_dir)
        if data.conda is not None:
            data.conda.write(file_dir)
        self._publish((UploadType.FILE, [runtime.to_file_model(file_dir)]))
        return runtime

    @backup()
    @traced()
    @async_io()
    def trace_log(self, data: LogData) -> list[BaseModel]:
        """
        追踪日志数据，日志数据比较特殊，一次可能好几行
        """
        logs = Log.from_log_data(data)
        self._publish((UploadType.LOG, [log.to_log_model() for log in logs]))
        return logs

    @once()
    @traced()
    def close_trace(self, success: bool, error: str = None, epoch: int = None):
        """
        停止日志跟踪，清理相关资源
        """
        # 上传错误日志
        if error is not None:
            log = Log.model_validate({"level": "ERROR", "message": error, "create_time": create_time(), "epoch": epoch})
            self._publish((UploadType.LOG, [log.to_log_model()]))
            # 备份日志
            self._f.write(log.to_record())
        # 停止worker
        if self._executor is not None:
            self._executor.shutdown(wait=True)
        # 停止上传线程池
        if self._pool is not None:
            self._pool.finish()
        # 写入结束标志
        footer = Footer.model_validate({"create_time": create_time(), "success": success})
        self._f.write(footer.to_record())
        self._f.ensure_flushed()
        self._f.close()

    def __new__(cls):
        """
        开启单例模式，确保同时只有一个 ProtoTransfer 实例存在，即代表一个实验会话
        创建会话之前，client 必须存在
        """
        if cls._instance is not None:
            return cls._instance
        # 检查是否已经创建 client
        try:
            get_client()
        except ValueError:
            raise ValueError("Client not initialized when creating ProtoTransfer instance")
        run_store = get_run_store()
        assert run_store.media_dir, "Media directory must be set before creating ProtoTransfer instance"
        assert run_store.file_dir, "File directory must be set before creating ProtoTransfer instance"
        assert run_store.backup_file, "Backup file must be set before creating ProtoTransfer instance"
        cls._run_store = run_store
        cls._instance = super(DataPorter, cls).__new__(cls)
        return cls._instance

    def open_for_sync(self):
        """ """
