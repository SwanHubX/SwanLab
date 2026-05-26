"""
@author: cunyue
@file: sync.py
@time: 2026/5/14 14:21
@description: sync 协议实现
"""

from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import DataStoreError
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
    DeliverSyncStartResponse,
)
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, ResumeMode, RunState, StartRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.api.experiment import get_experiment_summary, stop_experiment
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.metrics import RunMetrics
from swanlab.sdk.internal.core_python.pkg import builder, counter
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.core_python.transport import Transport
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.utils import (
    PrepareExperimentStartResult,
    generate_run_online_path,
    prepare_experiment_start,
)
from swanlab.sdk.internal.pkg import adapter, console, executor, safe
from swanlab.sdk.protocol.core import CoreSyncProtocol
from swanlab.sdk.typings.core_python.api.experiment import ResumeExperimentSummaryType

__all__ = ["CoreSyncPython"]


class CoreSyncPython(CoreSyncProtocol):
    def __init__(self):
        super().__init__()
        self._run_ctx: Optional[CoreContext] = None
        self._reader = DataStoreReader()
        self._read_executor = executor.EventLoopExecutor("SwanLab·Sync·Reader")
        self._transport: Optional[Transport] = None
        self._tracker: Optional[UploadTracker] = None
        self._start_request: Optional[DeliverSyncStartRequest] = None
        self._start_record: Optional[StartRecord] = None
        self._metrics: Optional[RunMetrics] = None
        # finish record 有可能不存在，此时Sync自己mock一下数据
        self._finish_record: Optional[FinishRecord] = None
        self._prepared_finish_record: Optional[FinishRecord] = None
        self._finish_records_queued = False
        self._finish_requested = False
        # 日志行数，仅大于此行数的才上传
        self._epoch = -1
        self._counter = counter.Counter(builder.DEC_NUM)

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        self._start_request = start_request
        self._ctx = CoreContext.from_proto(start_request.core_settings, mode="sync")
        # 安全地开启sync，DataStoreError 通常是run头文件损坏，属于已知错误
        with safe.block(message="Failed to open sync store with unexpected error"):
            try:
                self._reader.open(self._ctx.run_file)
                start_record_bytes = self._reader.scan()
                if start_record_bytes is None:
                    return DeliverSyncStartResponse(
                        success=False, message=f"Cannot scan records from {self._ctx.run_file}"
                    )
                # 第一个记录必然是 start_record 如果不是说明记录是损坏的，不符合生命周期设计，不过由于协议本身兜底，一般不会出现这种情况
                record = Record()
                record.ParseFromString(start_record_bytes)
                if not record.HasField("start"):
                    return DeliverSyncStartResponse(
                        success=False,
                        message=(
                            "Failed to sync this run because its log file is missing the required start record. "
                            "The file may be corrupted, incomplete, or not a valid SwanLab run file."
                        ),
                    )
                self._start_record = record.start
                return DeliverSyncStartResponse(success=True, message="success")
            except DataStoreError as e:
                return DeliverSyncStartResponse(success=False, message=str(e))
        return DeliverSyncStartResponse(success=False, message="unknown error")

    def deliver_sync_flush(self) -> DeliverSyncFlushResponse:
        assert self._start_request is not None and self._start_record is not None, (
            "Please run deliver_sync_start first before calling deliver_sync_flush"
        )

        result: Optional[PrepareExperimentStartResult] = None
        with safe.block(message="Failed to prepare sync experiment"):
            # 1. 通知后端，启动实验
            start_record = StartRecord()
            start_record.CopyFrom(self._start_record)
            if self._start_request.project:
                start_record.project = self._start_request.project
            if self._start_request.workspace:
                start_record.workspace = self._start_request.workspace
            if self._start_request.id:
                start_record.id = self._start_request.id
            # sync 实验统一开启宽松resume限制，resume设置为allow，后端自动创建新实验
            start_record.resume = ResumeMode.RESUME_MODE_ALLOW
            # sync的实验不同步实验颜色，否则会给用户一些奇怪的感觉：https://github.com/SwanHubX/SwanLab/issues/1434
            start_record.color = ""
            result = prepare_experiment_start(start_record)
            self._ctx.set_online_params(
                username=result.username,
                project=result.project,
                project_id=result.project_info["cuid"],
                experiment_id=result.experiment["cuid"],
            )

        if result is None:
            return DeliverSyncFlushResponse(success=False, message="Failed to prepare sync experiment")

        metrics_ready = False
        with safe.block(message="Failed to initialize sync metrics", write_to_tty=False):
            # 2. 如果是旧实验，向后端请求实验摘要
            summary: Optional[ResumeExperimentSummaryType] = None
            if not result.new_experiment:
                summary = get_experiment_summary(self._ctx.project_id, self._ctx.experiment_id)
            self._metrics, self._epoch, _, _ = RunMetrics.new(summary, ctx=self._ctx)
            metrics_ready = True

        if not metrics_ready:
            return DeliverSyncFlushResponse(success=False, message="Failed to initialize sync metrics")

        # 3. 开始读取本地文件并开始上传
        self._tracker = UploadTracker()
        self._tracker.set_state(CoreState.CORE_STATE_RUNNING)
        self._transport = Transport(self._ctx, tracker=self._tracker)
        self._read_executor.start(self._read_and_request_finish())
        self._transport.start()
        return DeliverSyncFlushResponse(success=True, message="success", path=generate_run_online_path(result))

    async def _read_and_request_finish(self):
        try:
            await self.read()
        finally:
            self._request_transport_finish()

    async def read(self):
        """
        异步读取本地文件
        """
        assert self._transport is not None, "Transport not set before reading records"
        assert self._metrics is not None, "Metrics not set before reading records"
        for record_bytes in self._reader:
            with safe.block(message="Failed to parse record"):
                record = Record()
                record.ParseFromString(record_bytes)
                # 记录开始和结束记录作为生命周期记录，不通过 _transport 上传
                if record.HasField("finish"):
                    if self._finish_record is None:
                        self._finish_record = record.finish
                    else:
                        console.error(
                            "Multiple finish records were found in the run file. SwanLab will use the first one."
                        )
                elif record.HasField("log"):
                    if record.log.epoch > self._epoch:
                        # 仅大于 console_epoch 的 record 才能上传
                        self._epoch = record.log.epoch
                        self._transport.put([record])
                    else:
                        console.debug(f"Log at epoch {record.log.epoch} was skipped syncing because it is too old.")
                elif record.HasField("column"):
                    # 如果 column 不存在于 _metrics 中则上传，并且定义此column
                    column = record.column
                    if self._metrics.get(column.column_key) is not None:
                        console.debug(f"Column '{column.column_key}' was skipped syncing because it already exists.")
                        continue
                    if column.column_type == ColumnType.COLUMN_TYPE_SCALAR:
                        self._metrics.define_scalar(key=column.column_key, column=column)
                    else:
                        self._metrics.define_media(
                            key=column.column_key, column=column, path=adapter.medium[column.column_type]
                        )
                    self._transport.put([record])
                elif record.HasField("scalar"):
                    # 如果scalar对应的column不存在于metrics（sync的实验本身就是resume的），则需要创建默认的column并上传
                    scalar = record.scalar
                    # 1. 创建列
                    metric = self._metrics.get(scalar.key)
                    if not metric:
                        column_record = builder.build_auto_column(self._ctx, scalar)
                        metric = self._metrics.define_scalar(key=scalar.key, column=column_record)
                        self._transport.put([builder.build_column_record(self._counter, column_record, positive=False)])
                    # 2. 检查类型是否匹配，判断指定的step是否允许写入
                    metric.ensure_type_match(scalar.type)
                    if metric.try_accept_step(scalar.step):
                        metric.update(scalar)
                        self._transport.put([record])
                    else:
                        console.debug(
                            f"Metric '{scalar.key}' at step {scalar.step} was skipped syncing because it is duplicate or too old."
                        )
                elif record.HasField("media"):
                    # 如果media对应的column不存在于metrics（sync的实验本身就是resume的），则需要创建默认的column并上传
                    media = record.media
                    # 1. 创建列
                    metric = self._metrics.get(media.key)
                    if not metric:
                        column_record = builder.build_auto_column(self._ctx, media)
                        metric = self._metrics.define_media(
                            key=media.key, column=column_record, path=self._ctx.media_dir / adapter.medium[media.type]
                        )
                        self._transport.put([builder.build_column_record(self._counter, column_record, positive=False)])
                    # 2. 检查类型是否匹配，判断指定的step是否允许写入
                    metric.ensure_type_match(media.type)
                    if metric.try_accept_step(media.step):
                        metric.update(media)
                        self._transport.put([record])
                    else:
                        console.debug(
                            f"Metric '{media.key}' at step {media.step} was skipped syncing because it is duplicate or too old."
                        )
                else:
                    # 其余的直接上传（覆盖式）
                    self._transport.put([record])

    def get_operation_stats(self) -> GetOperationStatsResponse:
        if self._tracker is None:
            return GetOperationStatsResponse(success=False, message="Transport not initialized yet")
        stats = self._tracker.snapshot()
        return GetOperationStatsResponse(success=True, message="success", stats=stats)

    def _get_finish_record(self) -> FinishRecord:
        if self._prepared_finish_record is not None:
            return self._prepared_finish_record

        finish_record = FinishRecord()
        if self._finish_record is not None:
            finish_record.CopyFrom(self._finish_record)
        else:
            # 找不到 Finish Record，一般情况为整个进程被Kill了，此时需要mock一个Finish Record
            finish_record.state = RunState.RUN_STATE_CRASHED
            ts = Timestamp()
            ts.GetCurrentTime()
            finish_record.finished_at.CopyFrom(ts)
            finish_record.error = (
                "Run process was interrupted or killed before writing a finish record. "
                "The run is marked as crashed during sync."
            )
        self._prepared_finish_record = finish_record
        return finish_record

    def _queue_finish_records(self) -> None:
        assert self._transport is not None, "Transport not set before queueing finish records"
        if self._finish_records_queued:
            return

        finish_record = self._get_finish_record()
        if finish_record.state != RunState.RUN_STATE_FINISHED:
            error_message = (
                finish_record.error
                if finish_record.error
                else "run failed with unknown error while finish_record.error is not set"
            )
            record = builder.build_log_record(
                self._counter,
                counter.Counter(self._epoch),
                LogRecord(
                    timestamp=finish_record.finished_at,
                    level=LogLevel.LOG_LEVEL_ERROR,
                    line=error_message,
                ),
                positive=False,
            )
            self._transport.put([record])
        self._finish_records_queued = True

    def _request_transport_finish(self) -> None:
        assert self._transport is not None, "Transport not set before requesting sync finish"
        if self._finish_requested:
            return

        self._queue_finish_records()
        self._transport.request_finish()
        self._finish_requested = True

    def confirm_sync_finish(self) -> ConfirmSyncFinishResponse:
        """
        确认同步完成，关闭文件流，上报最终实验状态
        """
        assert self._transport is not None, "Transport not set before confirming sync finish"
        # 1. 等待读任务完成，并构建结束记录
        self._read_executor.wait()
        self._read_executor.close()
        self._reader.close()
        finish_record = self._get_finish_record()
        # 3. 确保结束日志已入队，等待上传线程退出，然后上报最终实验状态
        try:
            self._request_transport_finish()
            self._transport.join(timeout=None)
            with safe.block(message="Failed to report run finish"):
                stop_experiment(
                    self._ctx.username,
                    self._ctx.project,
                    self._ctx.experiment_id,
                    state=finish_record.state,
                    finished_at=finish_record.finished_at,
                )
                self._pending_online_finish_record = None
                return ConfirmSyncFinishResponse(success=True, message="OK")
                # 如果仅仅是与后端同步出现问题，则换一个让用户安心一些的提示信息
            return ConfirmSyncFinishResponse(
                success=False, message="Failed to finish run, but all records have been uploaded to the server."
            )
        finally:
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)
