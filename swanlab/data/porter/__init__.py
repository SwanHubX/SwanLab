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
from typing import Optional, Literal, List, Union, Tuple

import wrapt

from swanlab.core_python import get_client
from swanlab.core_python.uploader import ColumnModel, ScalarModel, MediaModel, LogModel
from swanlab.core_python.uploader.thread import ThreadPool, UploadType
from swanlab.data.store import RunStore, get_run_store, reset_run_store
from swanlab.error import ValidationError
from swanlab.log.type import LogData
from swanlab.proto.v0 import Log, Header, Project, Experiment, Column, Metric, BaseModel, Runtime, Footer, Media, Scalar
from swanlab.toolkit import MetricInfo, ColumnInfo, RuntimeInfo, create_time, LogContent
from .datastore import DataStore
from .mounter import Mounter
from .utils import filter_metric, filter_epoch, filter_column

__all__ = ['DataPorter', 'Mounter']


def traced():
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        if instance is not None and getattr(instance, "_mode") != 1:
            raise RuntimeError("ProtoTransfer has not been started for tracing. Call oopen_for_trace() first.")
        return wrapped(*args, **kwargs)

    return wrapper


def synced():
    """
    标记需要同步上传的函数
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        if instance is not None and getattr(instance, "_mode") != 2:
            raise RuntimeError("ProtoTransfer has not been started for sync upload. Call open_for_sync() first.")
        return wrapped(*args, **kwargs)

    return wrapper


def backup():
    """
    标记需要备份的函数，并将函数的返回结果存储到备份文件中
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        result: Union[List[BaseModel], BaseModel] = wrapped(*args, **kwargs)
        if result is not None:
            try:
                f = getattr(instance, "_f")
                if isinstance(result, list):
                    [f.write(item.to_record()) for item in result]
                else:
                    f.write(result.to_record())
            except ValueError:
                # 写入备份文件失败，可能是因为没有开启备份模式或者备份文件未打开
                # TODO: 记录本地日志
                pass
            except TypeError:
                # 目前此错误会发生在与 transformers 集成时 ctrl + c，tqdm 没有立即退出
                # TODO: 记录本地日志
                pass
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
        # 当前模式
        # 0: 不开启任何模式
        # 1: 实验日志跟踪模式
        # 2: 同步上传模式，用于本地日志上报
        self._mode: Literal[0, 1, 2] = 0
        # 工作线程池，开启后一些方法会运行在子线程中
        self._executor: Optional[ThreadPoolExecutor] = None
        self._closed = False

        # ---------------------------------- 同步时用到的参数 ----------------------------------
        self._header: Optional[Header] = None
        self._project: Optional[Project] = None
        self._experiment: Optional[Experiment] = None
        self._logs: List[Log] = []
        self._runtime: Runtime = Runtime(
            conda_filename=None,
            requirements_filename=None,
            metadata_filename=None,
            config_filename=None,
        )
        self._columns: List[Column] = []
        self._scalars: List[Scalar] = []
        self._medias: List[Media] = []
        self._footer: Optional[Footer] = None

    def __new__(cls):
        """
        开启单例模式，确保同时只有一个 ProtoTransfer 实例存在，即代表一个实验会话
        创建会话之前，client 必须存在
        """
        if cls._instance is not None:
            raise RuntimeError("DataPorter instance already exists, cannot create a new one.")
        run_store = get_run_store()
        cls._run_store = run_store
        cls._instance = super(DataPorter, cls).__new__(cls)
        return cls._instance

    def _set_mode(self, mode: Literal[0, 1, 2]):
        """
        设置当前模式，只允许从0设置为其他模式
        """
        assert mode in [0, 1, 2], "Mode must be one of [0, 1, 2]."
        assert self._mode == 0, "DataPorter is already in use, cannot change mode."
        assert self._closed is False, "DataPorter has already rested, cannot change mode."
        self._mode = mode

    def open_for_trace(self, sync: bool = False, backend: Literal['go', 'python', 'none'] = 'python'):
        """
        开启日志传输器以实验日志跟踪，创建数据存储句柄
        此函数在同一实例中只能调用一次，后续调用会抛出异常
        :param sync: 是否同步上传数据，默认为 False，表示异步上传
        :param backend: 上传后端类型，默认为 'python'，可选 'go' 或 'none', none 表示不开启上传
        :raises RuntimeError: 如果已经开启了跟踪模式或者别的模式，则抛出异常
        """
        assert self._run_store.media_dir, "Media directory must be set before creating ProtoTransfer instance"
        assert self._run_store.file_dir, "File directory must be set before creating ProtoTransfer instance"
        assert self._run_store.backup_file, "Backup file must be set before creating ProtoTransfer instance"

        assert self._mode == 0, "DataPorter is already in use, cannot open for trace again."
        assert self._closed is False, "DataPorter has already rested, cannot open for trace again."

        self._f.open_for_write(self._run_store.backup_file)
        if backend == 'python':
            # 检查是否已经创建 client
            try:
                get_client()
            except ValueError:
                raise ValueError("Client not initialized when creating ProtoTransfer instance")
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
        # 项目
        project = self._run_store.project
        workspace = self._run_store.workspace
        visibility = self._run_store.visibility
        # 实验
        id = self._run_store.run_id
        name = self._run_store.run_name
        colors = self._run_store.run_colors
        description = self._run_store.description
        tags = self._run_store.tags
        self._f.write(
            Project.model_validate(
                {
                    "name": project,
                    "workspace": workspace,
                    "public": visibility,
                }
            ).to_record()
        )
        self._f.write(
            Experiment.model_validate(
                {
                    "id": id,
                    "name": name,
                    "colors": colors,
                    "description": description,
                    "tags": tags,
                }
            ).to_record()
        )
        self._set_mode(1)

    def _publish(self, *args, **kwargs):
        """
        发布数据到上传线程池，如果未开启线程池，则不执行任何操作
        """
        if self._pool is not None:
            self._pool.queue.put(*args, **kwargs)

    @async_io()
    @backup()
    @traced()
    def trace_column(self, data: ColumnInfo) -> BaseModel:
        """
        追踪列数据
        """
        column = Column.from_column_info(data)
        self._publish((UploadType.COLUMN, [column.to_column_model()]))
        return column

    @async_io()
    @backup()
    @traced()
    def trace_metric(self, data: MetricInfo) -> BaseModel:
        """
        追踪指标数据
        """
        assert data.error is None, "MetricInfo must not have error, if it has error, do not upload it."
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
                    if data.metric is None:
                        # TODO: 记录日志
                        raise ValueError("MetricInfo must have metric data when uploading media.")
                    # 写入数据
                    with open(os.path.join(path, data.metric["data"][i]), "wb") as f:
                        f.write(r.getvalue())
            self._publish((UploadType.MEDIA_METRIC, [media.to_media_model(data.swanlab_media_dir)]))
            return media

    @async_io()
    @backup()
    @traced()
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

    @async_io()
    @backup()
    @traced()
    def trace_log(self, data: LogData) -> List[BaseModel]:
        """
        追踪日志数据，日志数据比较特殊，一次可能好几行
        """
        logs = Log.from_log_data(data)
        self._publish((UploadType.LOG, [log.to_log_model() for log in logs]))
        return logs

    @traced()
    def close_trace(self, success: bool, error: str = None, epoch: int = None):
        """
        停止日志跟踪，清理相关资源
        """
        assert self._closed is False, "DataPorter has already rested, cannot close trace again."
        assert self._mode == 1, "DataPorter is not in trace mode, cannot close trace."

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
        self._closed = True
        DataPorter._reset()

    @classmethod
    def _reset(cls):
        """重置单例实例，允许创建新的DataPorter"""
        cls._instance = None
        cls._run_store = None

    def open_for_sync(self, run_dir: str, backend: Literal['python', 'go'] = 'python') -> "DataPorter":
        """
        开启同步模式，此函数应该与 __enter__() 一起使用
        """
        self._run_store.run_dir = run_dir
        assert self._run_store.media_dir, "Media directory must be set before creating ProtoTransfer instance"
        assert self._run_store.file_dir, "File directory must be set before creating ProtoTransfer instance"
        assert self._run_store.backup_file, "Backup file must be set before creating ProtoTransfer instance"
        backup_file = self._run_store.backup_file
        assert os.path.isfile(backup_file), f"Backup file {backup_file} does not exist."
        self._f.open_for_scan(backup_file)
        if backend == 'python':
            self._pool = ThreadPool()
        elif backend == 'go':
            raise NotImplementedError("swanlab-core is not ready yet.")
        else:
            raise ValueError(f"Unsupported backend for sync: {backend}")
        self._set_mode(2)
        return self

    @synced()
    def __enter__(self):
        """
        开启本地数据备份和上传线程池
        使用方式为：
        >>> backup_file = "path/to/backup/file"
        >>> with DataPorter().open_for_sync(backup_file) as porter:
        >>>    porter.parse()
        >>>    ... # 创建实验
        >>>    porter.synchronize()
        """
        assert self._mode == 2, "DataPorter is already in use, cannot open for sync."
        assert self._closed is False, "DataPorter has already rested, cannot open for sync."
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        reset_run_store()
        self._closed = True
        DataPorter._reset()
        # 出现错误则抛出异常
        if exc_type is not None:
            raise exc_val

    @synced()
    def parse(self) -> Tuple[Project, Experiment]:
        """
        解析备份文件中的记录，必须在 open_for_sync() 后调用
        """
        try:
            for record in self._f:
                if record is not None:
                    self._parse_record(record)
        except ValidationError:
            # 略过坏的记录，直接退出
            # TODO 写入日志文件
            pass
        # 检查是否所有必要的记录都已解析
        assert self._header is not None, "Header not parsed"
        # 检查备份文件
        assert (
            self._header.backup_type == "DEFAULT"
        ), f"Backup type mismatch: {self._header.backup_type}, please update your swanlab package."
        assert self._project is not None, "Project not parsed"
        assert self._experiment is not None, "Experiment not parsed"
        return self._project, self._experiment

    def _parse_record(self, data: str):
        assert self._mode == 2, "Must parse records in sync mode (mode=2)."
        assert self._closed is False, "DataPorter has already rested, cannot parse."
        record = BaseModel.from_record(data)
        if isinstance(record, Header):
            assert self._header is None, "Header already parsed"
            self._header = record
            return
        if isinstance(record, Project):
            assert self._project is None, "Project already parsed"
            self._project = record
            return
        if isinstance(record, Experiment):
            assert self._experiment is None, "Experiment already parsed"
            self._experiment = record
            return
        if isinstance(record, Log):
            self._logs.append(record)
            return
        if isinstance(record, Runtime):
            if record.conda_filename is not None:
                self._runtime.conda_filename = record.conda_filename
            if record.requirements_filename is not None:
                self._runtime.requirements_filename = record.requirements_filename
            if record.metadata_filename is not None:
                self._runtime.metadata_filename = record.metadata_filename
            if record.config_filename is not None:
                self._runtime.config_filename = record.config_filename
            return
        if isinstance(record, Column):
            self._columns.append(record)
            return
        if isinstance(record, Scalar):
            self._scalars.append(record)
            return
        if isinstance(record, Media):
            self._medias.append(record)
            return
        if isinstance(record, Footer):
            assert self._footer is None, "Footer already parsed"
            self._footer = record
            return
        raise ValueError(f"Unknown record type: {type(record)}")

    def close(self):
        """
        关闭实例，此函数用于没有开启 trace 和 sync 时关闭实例
        """
        assert self._closed is False, "DataPorter has already rested, cannot close."
        assert self._mode == 0, "DataPorter is in use, cannot simply close."
        DataPorter._reset()

    @synced()
    def synchronize(self):
        """
        同步上传数据到 SwanLab 服务器，必须在 open_for_sync() 后调用
        NOTE: 执行此函数后，实验会话被关闭
        """
        assert self._mode == 2, "Must synchronize in sync mode (mode=2)."
        assert self._closed is False, "DataPorter has already rested, cannot synchronize."
        # 同步上传数据到 SwanLab 服务器
        # 1. 上传文件（配置、运行时）
        self._publish((UploadType.FILE, [self._runtime.to_file_model(file_dir=self._run_store.file_dir)]))
        # 2. 上传列信息
        columns = filter(self._filter_column_by_key, [column.to_column_model() for column in self._columns])
        self._publish((UploadType.COLUMN, list(columns)))
        # 3. 上传标量指标
        scalars = filter(self._filter_scalar_by_step, [scalar.to_scalar_model() for scalar in self._scalars])
        self._publish((UploadType.SCALAR_METRIC, list(scalars)))
        # 4. 上传媒体指标
        medias = filter(
            self._filter_media_by_step, [media.to_media_model(self._run_store.media_dir) for media in self._medias]
        )
        self._publish((UploadType.MEDIA_METRIC, list(medias)))
        # 5. 上传日志
        logs = [log.to_log_model() for log in self._logs]
        for ls in logs:
            ls['contents'] = list(filter(self._filter_log_by_epoch, ls['contents']))
        self._publish((UploadType.LOG, logs))
        # 6. 等待上传完毕
        self._pool.finish()
        # 7. 同步实验状态
        client = get_client()
        if self._footer is None:
            client.update_state(success=False)
        else:
            client.update_state(self._footer.success, self._footer.create_time)

    def _filter_column_by_key(self, column: ColumnModel) -> bool:
        """
        筛选列指标，排除已经上传的列，此函数用于 filter 高阶函数
        :param column: 列指标数据
        :return 是否需要被保留
        """
        return filter_column(column.key, self._run_store.metrics)

    def _filter_scalar_by_step(self, metric: ScalarModel) -> bool:
        """
        筛选标量指标数据，排除已经上传的指标
        :param metric: 指标数据
        :return: 是否需要被保留
        """
        return filter_metric(metric.key, metric.step, self._run_store.metrics)

    def _filter_media_by_step(self, metric: MediaModel) -> bool:
        """
        筛选媒体数据，排除已经上传的媒体
        :param metric : 媒体数据
        :return: 是否需要被保留
        """
        return filter_metric(metric.key, metric.step, self._run_store.metrics)

    def _filter_log_by_epoch(self, log: LogContent) -> LogModel:
        """
        筛选日志数据，排除已经上传的日志
        :param log: 日志数据
        :return: 是否需要被保留
        """
        return filter_epoch(log['epoch'], self._run_store.log_epoch)
