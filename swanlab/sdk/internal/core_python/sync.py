"""
@author: cunyue
@file: sync.py
@time: 2026/5/14 14:21
@description: sync 协议实现
"""

from typing import Optional, Tuple, Union

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
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, ResumeMode, RunState, StartRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.api.experiment import get_experiment_summary, stop_experiment
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.metrics import RunMetrics
from swanlab.sdk.internal.core_python.pkg import builder, counter
from swanlab.sdk.internal.core_python.pkg.range import SyncRangeFilter
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
        # finish record 只来自本地文件，缺失表示运行未正常结束或仍在写
        self._finish_record: Optional[FinishRecord] = None
        self._finish_records_queued = False
        self._finish_requested = False
        # 日志行数，仅大于此行数的才上传
        self._epoch = -1
        self._counter = counter.Counter(builder.DEC_NUM)
        # 范围同步过滤器，默认关闭表示全量同步
        self._filter = SyncRangeFilter()
        # 本次 sync 是否新建了实验，用于决定缺失 finish record 时是否按 CRASHED 收尾
        self._new_experiment = False

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        self._start_request = start_request
        self._filter = SyncRangeFilter.from_proto(start_request.range)
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
        # 记录服务端拒绝原因（403、COPY 实验等），随失败响应返回给调用方
        prepare_errors: list[BaseException] = []
        with safe.block(message="Failed to prepare sync experiment", on_error=prepare_errors.append):
            # 1. 通知后端，创建或恢复实验
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
            # 已有实验会通过 Server 的 resume 流程置回 RUNNING，
            # 使本次回填写入落在有效查询窗口内（结束后再按真实 FinishRecord 收尾）
            result = prepare_experiment_start(start_record)
            self._ctx.set_online_params(
                username=result.username,
                project=result.project,
                project_id=result.project_data["cuid"],
                project_version=result.project_data.get("version", None),
                experiment_id=result.experiment["cuid"],
            )

        if result is None:
            detail = f": {prepare_errors[0]}" if prepare_errors else ""
            return DeliverSyncFlushResponse(success=False, message=f"Failed to prepare sync experiment{detail}")

        self._new_experiment = result.new_experiment

        metrics_ready = False
        with safe.block(message="Failed to initialize sync metrics", write_to_tty=False):
            # 2. 如果是旧实验，向后端请求实验摘要
            summary: Optional[ResumeExperimentSummaryType] = None
            if not result.new_experiment:
                summary = get_experiment_summary(
                    self._ctx.project_id,
                    self._ctx.experiment_id,
                    created_at=result.experiment["createdAt"],
                )
            self._metrics, self._epoch, _, _ = RunMetrics.new(summary, ctx=self._ctx)
            if self._filter.enabled:
                # 范围同步用于回放/恢复历史数据，必须允许写回已存在的 step
                self._metrics.allow_replay()
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
                self._handle_record(record)

    def _handle_record(self, record: Record) -> None:
        """按记录类型分发处理。"""
        assert self._transport is not None, "Transport not set before reading records"
        # 记录开始和结束记录作为生命周期记录，不通过 _transport 上传
        if record.HasField("finish"):
            self._handle_finish(record.finish)
        elif record.HasField("log"):
            self._handle_log(record)
        elif record.HasField("column"):
            self._handle_column(record)
        elif record.HasField("scalar"):
            self._handle_metric(record, record.scalar, is_media=False)
        elif record.HasField("media"):
            self._handle_metric(record, record.media, is_media=True)
        else:
            # 其余的直接上传（覆盖式）
            self._transport.put([record])

    def _handle_finish(self, finish: FinishRecord) -> None:
        if self._finish_record is None:
            self._finish_record = finish
        else:
            console.error("Multiple finish records were found in the run file. SwanLab will use the first one.")

    def _handle_log(self, record: Record) -> None:
        assert self._transport is not None, "Transport not set before reading records"
        log = record.log
        if self._filter.enabled:
            # step 范围只约束指标；终端日志仅按时间范围回放
            if not self._filter.has_time:
                return
            if not self._filter.match(timestamp=self._pick_timestamp(log.timestamp, record.timestamp)):
                return
            self._epoch = max(self._epoch, log.epoch)
            self._transport.put([record])
        elif log.epoch > self._epoch:
            # 仅大于 console_epoch 的 record 才能上传
            self._epoch = log.epoch
            self._transport.put([record])
        else:
            console.debug(f"Log at epoch {log.epoch} was skipped syncing because it is too old.")

    def _handle_column(self, record: Record) -> None:
        assert self._transport is not None, "Transport not set before reading records"
        assert self._metrics is not None, "Metrics not set before reading records"
        # 如果 column 不存在于 _metrics 中则上传，并且定义此column
        column = record.column
        if self._metrics.get(column.column_key) is not None:
            console.debug(f"Column '{column.column_key}' was skipped syncing because it already exists.")
            return
        if column.column_type == ColumnType.COLUMN_TYPE_SCALAR:
            self._metrics.define_scalar(key=column.column_key, column=column)
        else:
            self._metrics.define_media(key=column.column_key, column=column, path=adapter.medium[column.column_type])
        self._transport.put([record])

    def _handle_metric(self, record: Record, data: Union[ScalarRecord, MediaRecord], *, is_media: bool) -> None:
        """处理 scalar/media 数据记录，按需自动补建列。"""
        assert self._transport is not None, "Transport not set before reading records"
        assert self._metrics is not None, "Metrics not set before reading records"
        if self._filter.enabled and not self._filter.match(
            step=data.step, timestamp=self._pick_timestamp(data.timestamp, record.timestamp)
        ):
            return
        # 1. 创建列
        metric = self._metrics.get(data.key)
        if not metric:
            column = builder.build_auto_column(self._ctx, data)
            if is_media:
                metric = self._metrics.define_media(
                    key=data.key, column=column, path=self._ctx.media_dir / adapter.medium[data.type]
                )
            else:
                metric = self._metrics.define_scalar(key=data.key, column=column)
            self._transport.put([builder.build_column_record(self._counter, column, positive=False)])
        # 2. 检查类型是否匹配，判断指定的step是否允许写入
        metric.ensure_type_match(data.type)
        if metric.try_accept_step(data.step):
            metric.update(data)
            self._transport.put([record])
        else:
            console.debug(
                f"Metric '{data.key}' at step {data.step} was skipped syncing because it is duplicate or too old."
            )

    @staticmethod
    def _pick_timestamp(payload: Timestamp, record: Timestamp) -> Timestamp:
        """优先使用业务时间（指标/日志时间），未设置时回退到记录落盘时间。"""
        return payload if payload.seconds or payload.nanos else record

    def get_operation_stats(self) -> GetOperationStatsResponse:
        if self._tracker is None:
            return GetOperationStatsResponse(success=False, message="Transport not initialized yet")
        stats = self._tracker.snapshot()
        return GetOperationStatsResponse(success=True, message="success", stats=stats)

    def _queue_finish_records(self) -> None:
        assert self._transport is not None, "Transport not set before queueing finish records"
        if self._finish_records_queued:
            return

        # 范围同步只回放选定时间片，不注入原实验结束时的错误日志
        if self._filter.enabled:
            self._finish_records_queued = True
            return

        finish_record = self._finish_record
        # 没有真实 FinishRecord 时不推断运行失败，也不生成“进程中断”日志
        if finish_record is None:
            self._finish_records_queued = True
            return

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

    def _resolve_finish_state(self) -> Optional[Tuple["RunState", Timestamp]]:
        """确定本次 sync 是否上报实验终态。

        - 有真实 FinishRecord：按记录中的状态收尾；
        - 没有记录且实验是本次新建：按 CRASHED 收尾（离线补传的既有行为）；
        - 没有记录且恢复到已有实验：只结束本次同步，实验保持 RUNNING，
          避免误停其他终端仍在训练的实验。
        """
        if self._finish_record is not None:
            return self._finish_record.state, self._finish_record.finished_at
        if self._new_experiment:
            finished_at = Timestamp()
            finished_at.GetCurrentTime()
            return RunState.RUN_STATE_CRASHED, finished_at
        return None

    def confirm_sync_finish(self) -> ConfirmSyncFinishResponse:
        """
        确认同步完成，关闭文件流，上报最终实验状态
        """
        assert self._transport is not None, "Transport not set before confirming sync finish"
        try:
            # 1. 等待读任务完成
            self._read_executor.wait()
            finish_state = self._resolve_finish_state()
            # 2. 触发结束日志入队，等待上传线程排空并上报实验状态
            self._request_transport_finish()
            # join 只保证线程退出，finish 才能确认缓冲区与未确认记录已排空
            if not self._transport.finish(timeout=None):
                return ConfirmSyncFinishResponse(success=False, message=self._transport.finish_summary())
            if finish_state is None:
                return ConfirmSyncFinishResponse(success=True, message="OK")
            state, finished_at = finish_state
            with safe.block(message="Failed to report run finish"):
                stop_experiment(
                    self._ctx.username,
                    self._ctx.project,
                    self._ctx.experiment_id,
                    state=state,
                    finished_at=finished_at,
                )
                return ConfirmSyncFinishResponse(success=True, message="OK")
                # 如果仅仅是与后端同步出现问题，则换一个让用户安心一些的提示信息
            return ConfirmSyncFinishResponse(
                success=False, message="Failed to finish run, but all records have been uploaded to the server."
            )
        finally:
            self._read_executor.close()
            self._reader.close()
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)
