"""
@author: cunyue
@file: core.py
@time: 2026/5/14 14:21
@description: core 协议实现，SwanLab Core Python 版本，封装SwanLab云端版核心业务，包括：
1. 提供http客户端，用于与SwanLab云端API进行交互。
2. 提供rpc封装函数，以rpc方式调用SwanLab云端API。
3. 提供上传线程，在另一个线程执行上传任务。
4. 存储指标上下文，便于分布式训练时的指标同步。
...

实现 CoreProtocol，当前为纯 Python 实现。
未来由 swanlab-core（Go 二进制）替代时，此模块整体被替换，
BackgroundConsumer 等调用方无需修改。


Core 同时需要根据不同模式处理不同的业务，这是设计模式决定的
值得说明的是，在当前的上层设计中，upsert 方法在 disabled 模式下永远不会触发，但是考虑到设计完整性，我们增加了相关业务逻辑判断
"""

from typing import List, Optional

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.env.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, RunState, StartRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.api.experiment import (
    get_experiment_summary,
    stop_experiment,
)
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.heartbeat import Heartbeat
from swanlab.sdk.internal.core_python.metrics import RunMetrics
from swanlab.sdk.internal.core_python.pkg import builder, counter
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.transport import Transport
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.utils import generate_run_online_path, prepare_experiment_start
from swanlab.sdk.internal.core_python.watcher import FileWatcher, create_save_links
from swanlab.sdk.internal.pkg import adapter, console, safe
from swanlab.sdk.protocol import CoreProtocol
from swanlab.sdk.typings.core_python.api.experiment import ResumeExperimentSummaryType
from swanlab.sdk.typings.run import ModeType

__all__ = ["CorePython"]


class CorePython(CoreProtocol):
    """
    CoreProtocol 的 Python 实现。
    由 Run 在初始化时构造并注入给 BackgroundConsumer
    """

    def __init__(self, mode: ModeType):
        super().__init__(mode)
        self._started: bool = False
        self._run_ctx: Optional[CoreContext] = None
        self._store: Optional[DataStoreWriter] = None
        self._transport: Optional[Transport] = None
        # record 构建计数器
        self._counter = counter.Counter()
        # console 行计数器
        self._epoch = counter.Counter()
        # 指标上下文管理器
        self._metrics: Optional[RunMetrics] = None
        # 在线模式心跳，定时向后端发送心跳以维持在线状态
        self._heartbeat: Optional[Heartbeat] = None
        # 在线模式上传跟踪器，记录上传队列状态，供进度展示使用
        self._tracker: Optional[UploadTracker] = None
        # finish 时暂存的记录，等待 confirm_run_finish 时上报，用于online模式的两阶段 finish 设计
        self._pending_online_finish_record: Optional[FinishRecord] = None
        # save 相关
        self._watcher = FileWatcher(on_change=self._on_file_changed)
        self._pending_end_saves: List[Record] = []

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    # ---------------------------------- 实验开始 ----------------------------------

    def deliver_run_start(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        if self._started:
            raise RuntimeError("Failed to start run: already started")
        resp = super().deliver_run_start(start_request)
        self._started = resp.success
        return resp

    def _start_store(self, resp: DeliverRunStartResponse):
        self._store = DataStoreWriter()
        self._store.open(str(self._ctx.run_file))
        record = builder.build_start_record(resp.run)
        self._store.write(record.SerializeToString())

    def _start_without_online(self, start_request: DeliverRunStartRequest, message: str) -> DeliverRunStartResponse:
        self._ctx = CoreContext.from_proto(start_request.core_settings)
        self._metrics, console_epoch, global_step, global_system_step = RunMetrics.new(None, ctx=self._ctx)
        resp = DeliverRunStartResponse(
            success=True,
            message=message,
            run=start_request.start_record,
            new_experiment=True,
            global_system_step=global_system_step,
            global_step=global_step,
        )
        self._epoch.reset(console_epoch)
        self._start_store(resp)
        return resp

    def _start_when_local(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._start_without_online(start_request, "OK, but use local")

    def _start_when_offline(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._start_without_online(start_request, "OK, but use offline")

    def _start_when_online(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        self._ctx = CoreContext.from_proto(start_request.core_settings)
        resp = self._report_run_start(start_request.start_record)
        self._start_store(resp)
        self._tracker = UploadTracker()
        self._tracker.set_state(CoreState.CORE_STATE_RUNNING)
        self._transport = Transport(ctx=self._ctx, tracker=self._tracker)
        self._heartbeat = Heartbeat(self._ctx.experiment_id)
        self._heartbeat.start()
        return resp

    def _report_run_start(self, record: StartRecord) -> DeliverRunStartResponse:
        """
        运行开始
        :param record: 运行开始记录
        :return: 运行开始响应
        """
        run_info = prepare_experiment_start(record)
        # 3. 记录必要字段
        self._ctx.set_online_params(
            username=run_info.username,
            project=run_info.project,
            project_id=run_info.project_info["cuid"],
            experiment_id=run_info.experiment["cuid"],
        )
        # 3. resume 时，向后端获取数据
        summary: Optional[ResumeExperimentSummaryType] = None
        if not run_info.new_experiment:
            summary = get_experiment_summary(self._ctx.project_id, self._ctx.experiment_id)
        self._metrics, console_epoch, global_step, global_system_step = RunMetrics.new(summary, ctx=self._ctx)
        self._epoch.reset(console_epoch)
        # 4. 构建记录
        start_record = StartRecord()
        start_record.CopyFrom(record)
        start_record.name = run_info.name
        start_record.color = run_info.color
        start_record.resume = record.resume
        start_record.project = run_info.project
        start_record.workspace = run_info.username
        return DeliverRunStartResponse(
            success=True,
            message="OK",
            path=generate_run_online_path(run_info),
            run=start_record,
            name=run_info.experiment.get("name"),
            global_step=global_step,
            global_system_step=global_system_step,
            new_experiment=run_info.new_experiment,
        )

    # ---------------------------------- 数据上报 ----------------------------------

    def _store_records(self, records: List[Record]) -> None:
        """将一组 Record 写入本地存储"""
        assert self._store is not None, "store must be initialized before upsert"
        for record in records:
            self._store.write(record.SerializeToString())

    def _transport_put(self, records: List[Record]) -> None:
        """将一组 Record 推送到上传队列"""
        assert self._transport is not None, "transport must be initialized before upsert"
        self._transport.put(records)

    # ---- upsert_columns ----

    def _upsert_columns_to_metrics(self, columns: List[ColumnRecord]):
        assert self._metrics is not None, "metrics must be initialized before upsert columns"
        records: List[Record] = []
        for column in columns:
            # 1. 已存在则跳过
            if self._metrics.get(column.column_key):
                handle = console.debug if column.column_class == ColumnClass.COLUMN_CLASS_SYSTEM else console.warning
                handle(f"Column {column.column_key} already has been defined, skipping")
                continue
            # 2. 否则定义指标，根据类型不同，定义不同的指标
            if column.column_type == ColumnType.COLUMN_TYPE_SCALAR:
                self._metrics.define_scalar(key=column.column_key, column=column)
            else:
                self._metrics.define_media(
                    key=column.column_key, column=column, path=adapter.medium[column.column_type]
                )
            records.append(builder.build_column_record(self._counter, column))
        self._store_records(records)
        return records

    def _upsert_columns_when_local(self, columns: List[ColumnRecord]) -> None:
        self._upsert_columns_to_metrics(columns)

    def _upsert_columns_when_offline(self, columns: List[ColumnRecord]) -> None:
        self._upsert_columns_to_metrics(columns)

    def _upsert_columns_when_online(self, columns: List[ColumnRecord]) -> None:
        records = self._upsert_columns_to_metrics(columns)
        self._transport_put(records)

    # ---- upsert_scalars ----

    def _upsert_scalars_to_metrics(self, scalars_list: List[ScalarRecord]) -> List[Record]:
        """
        将一组 ScalarRecord 转换为一组 Record，并更新指标上下文
        """
        assert self._metrics is not None, "metrics must be initialized before upsert scalars"
        records: List[Record] = []
        for scalar in scalars_list:
            with safe.block(message="Failed to upsert scalar metric ctx"):
                # 1. 如果指标未定义，则定义此指标
                metric = self._metrics.get(scalar.key)
                if not metric:
                    column_record = builder.build_auto_column(self._ctx, scalar)
                    metric = self._metrics.define_scalar(key=scalar.key, column=column_record)
                    records.append(builder.build_column_record(self._counter, column_record))
                # 2. 检查类型是否匹配，判断指定的step是否允许写入
                metric.ensure_type_match(scalar.type)
                if metric.try_accept_step(scalar.step):
                    metric.update(scalar)
                    records.append(builder.build_scalar_record(self._counter, scalar))
                else:
                    console.debug(
                        f"Metric '{scalar.key}' at step {scalar.step} was skipped because it is duplicate or too old."
                    )
        # 持久化存储record
        self._store_records(records)
        return records

    def _upsert_scalars_when_local(self, scalars: List[ScalarRecord]) -> None:
        self._upsert_scalars_to_metrics(scalars)

    def _upsert_scalars_when_offline(self, scalars: List[ScalarRecord]) -> None:
        self._upsert_scalars_to_metrics(scalars)

    def _upsert_scalars_when_online(self, scalars: List[ScalarRecord]) -> None:
        records = self._upsert_scalars_to_metrics(scalars)
        self._transport_put(records)

    # ---- upsert_media -----

    def _upsert_media_to_metrics(self, media_list: List[MediaRecord]) -> List[Record]:
        assert self._metrics is not None, "metrics must be initialized before upsert media"
        records: List[Record] = []
        for media in media_list:
            with safe.block(message="Failed to upsert media metric ctx"):
                # 1. 如果指标未定义，则定义此指标
                metric = self._metrics.get(media.key)
                if not metric:
                    column_record = builder.build_auto_column(self._ctx, media)
                    metric = self._metrics.define_media(
                        key=media.key, column=column_record, path=self._ctx.media_dir / adapter.medium[media.type]
                    )
                    records.append(builder.build_column_record(self._counter, column_record))
                # 2. 检查类型是否匹配，判断指定的step是否允许写入
                metric.ensure_type_match(media.type)
                if metric.try_accept_step(media.step):
                    metric.update(media)
                    records.append(builder.build_media_record(self._counter, media))
                else:
                    console.debug(
                        f"Metric '{media.key}' at step {media.step} was skipped because it is duplicate or too old."
                    )
            # 持久化存储record
        self._store_records(records)
        return records

    def _upsert_media_when_local(self, media: List[MediaRecord]) -> None:
        self._upsert_media_to_metrics(media)

    def _upsert_media_when_offline(self, media: List[MediaRecord]) -> None:
        self._upsert_media_to_metrics(media)

    def _upsert_media_when_online(self, media: List[MediaRecord]) -> None:
        records = self._upsert_media_to_metrics(media)
        self._transport_put(records)

    # ---- upsert_logs ----

    def _upsert_logs_when_local(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)

    def _upsert_logs_when_offline(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)

    def _upsert_logs_when_online(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)
        self._transport_put(records)

    # ---- upsert_configs ----

    def _upsert_configs_when_local(self, configs: List[ConfigRecord]) -> None:
        records = [builder.build_config_record(c) for c in configs]
        self._store_records(records)

    def _upsert_configs_when_offline(self, configs: List[ConfigRecord]) -> None:
        records = [builder.build_config_record(c) for c in configs]
        self._store_records(records)

    def _upsert_configs_when_online(self, configs: List[ConfigRecord]) -> None:
        records = [builder.build_config_record(c) for c in configs]
        self._store_records(records)
        self._transport_put(records)

    # ---- upsert_requirements ----

    def _upsert_requirements_when_local(self, requirements: List[RequirementsRecord]) -> None:
        records = [builder.build_requirements_record(r.timestamp) for r in requirements]
        self._store_records(records)

    def _upsert_requirements_when_offline(self, requirements: List[RequirementsRecord]) -> None:
        records = [builder.build_requirements_record(r.timestamp) for r in requirements]
        self._store_records(records)

    def _upsert_requirements_when_online(self, requirements: List[RequirementsRecord]) -> None:
        records = [builder.build_requirements_record(r.timestamp) for r in requirements]
        self._store_records(records)
        self._transport_put(records)

    # ---- upsert_conda ----

    def _upsert_conda_when_local(self, conda: List[CondaRecord]) -> None:
        records = [builder.build_conda_record(c.timestamp) for c in conda]
        self._store_records(records)

    def _upsert_conda_when_offline(self, conda: List[CondaRecord]) -> None:
        records = [builder.build_conda_record(c.timestamp) for c in conda]
        self._store_records(records)

    def _upsert_conda_when_online(self, conda: List[CondaRecord]) -> None:
        records = [builder.build_conda_record(c.timestamp) for c in conda]
        self._store_records(records)
        self._transport_put(records)

    # ---- upsert_metadata ----

    def _upsert_metadata_when_local(self, metadata: List[MetadataRecord]) -> None:
        records = [builder.build_metadata_record(m.timestamp) for m in metadata]
        self._store_records(records)

    def _upsert_metadata_when_offline(self, metadata: List[MetadataRecord]) -> None:
        records = [builder.build_metadata_record(m.timestamp) for m in metadata]
        self._store_records(records)

    def _upsert_metadata_when_online(self, metadata: List[MetadataRecord]) -> None:
        records = [builder.build_metadata_record(m.timestamp) for m in metadata]
        self._store_records(records)
        self._transport_put(records)

    # ---- upsert_saves ----
    def _create_save_with_links(self, saves: List[SaveRecord]) -> None:
        linked = create_save_links(saves, self._ctx.files_dir)
        if linked > 0:
            console.info(
                f"Symlinked {linked} files into the SwanLab run directory; call swanlab.save again to sync new files."
            )

    def _upsert_saves_when_local(self, saves: List[SaveRecord]) -> None:
        self._create_save_with_links(saves)
        records = [builder.build_save_record(self._counter, s) for s in saves]
        self._store_records(records)
        self._watcher.register_live_watches(saves, self._ctx.files_dir)

    def _upsert_saves_when_offline(self, saves: List[SaveRecord]) -> None:
        self._create_save_with_links(saves)
        records = [builder.build_save_record(self._counter, s) for s in saves]
        self._store_records(records)
        self._watcher.register_live_watches(saves, self._ctx.files_dir)

    def _upsert_saves_when_online(self, saves: List[SaveRecord]) -> None:
        self._create_save_with_links(saves)
        records = [builder.build_save_record(self._counter, s) for s in saves]
        self._store_records(records)
        self._watcher.register_live_watches(saves, self._ctx.files_dir)
        # "end" policy saves: store locally but defer cloud upload until finish
        transport_records: List[Record] = []
        for save, record in zip(saves, records):
            if save.policy == adapter.policy["end"]:
                self._pending_end_saves.append(record)
            else:
                transport_records.append(record)
        if transport_records:
            self._transport_put(transport_records)

    # ---- file watcher ----

    def _on_file_changed(self, save_record: SaveRecord) -> None:
        """FileWatcher 回调：持久化 + 可选上传。"""
        if not self._started or self._store is None:
            return
        records = [builder.build_save_record(self._counter, save_record)]
        self._store_records(records)
        if self._mode == "online" and self._transport is not None:
            self._transport_put(records)

    # ---------------------------------- fork 方法 ----------------------------------

    def fork(self) -> "CorePython":
        raise RuntimeError("CorePython.fork() should not be called (designed for swanlab-core).")

    # ---------------------------------- finish 方法 ----------------------------------

    def deliver_run_finish(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        if not self._started:
            return DeliverRunFinishResponse(success=False, message="Failed to finish run: not started")
        resp = super().deliver_run_finish(finish_request)
        self._started = False
        return resp

    def _store_finish(self, finish_record: FinishRecord) -> Optional[Record]:
        assert self._store is not None, "store must be initialized before shutdown"
        # 1. 构建结束记录并写入存储
        record = builder.build_finish_record(finish_record)
        record.finish.CopyFrom(finish_record)
        log_record: Optional[Record] = None
        if finish_record.state != RunState.RUN_STATE_FINISHED:
            error_message = (
                finish_record.error
                if finish_record.error
                else "run failed with unknown error while finish_record.error is not set"
            )
            log_record = builder.build_log_record(
                self._counter,
                self._epoch,
                LogRecord(
                    timestamp=finish_record.finished_at,
                    level=LogLevel.LOG_LEVEL_ERROR,
                    line=error_message,
                ),
            )
            # 不将 log_record 写入 store 中，一方面具体的报错信息存储在 finish_record 中
            # 另一方面因为这个 record 也是为了适应后端“报错信息写在CH”的设计
        # 2. 关闭存储、文件监视器等本地资源，停止接受新的记录
        self._store.write(record.SerializeToString())
        self._store.close()
        self._store = None
        self._watcher.stop()
        return log_record

    def _finish_when_local(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        self._store_finish(finish_request.finish_record)
        return DeliverRunFinishResponse(success=True, message="OK, but use local")

    def _finish_when_offline(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        self._store_finish(finish_request.finish_record)
        return DeliverRunFinishResponse(success=True, message="OK, but use offline")

    def _finish_when_online(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        """Online finish 分两阶段完成：
        1. deliver_run_finish：本地持久化 + 暂存 finish_record + 通知 transport 排空。
        2. confirm_run_finish（由 Run.finish 在进度展示后调用）：等待 transport 排空完成，再向
           后端上报最终实验状态。两阶段设计允许中间插入进度轮询展示。
        """
        assert self._transport is not None, "transport must be initialized before finishing"
        record = self._store_finish(finish_request.finish_record)
        # 2. 发送 error log 和暂存的 end-policy save 记录
        if record is not None:
            self._transport_put([record])
        if self._pending_end_saves:
            self._transport_put(self._pending_end_saves)
            self._pending_end_saves.clear()
        # 3. 暂存 finish_record，等待 confirm_run_finish 做最终确认
        self._pending_online_finish_record = finish_request.finish_record
        # 4. 通知 Transport 开始排空（非阻塞）
        self._transport.request_finish()
        return DeliverRunFinishResponse(success=True, message="OK")

    # ---------------------------------- 进度查询 ----------------------------------

    def _get_operation_when_oline(self) -> GetOperationStatsResponse:
        assert self._tracker is not None, "tracker must be initialized before get_operation_stats"
        stats = self._tracker.snapshot()
        return GetOperationStatsResponse(
            success=True,
            message="OK",
            stats=stats,
        )

    # ---------------------------------- 确认完成 ----------------------------------

    def _confirm_finish_when_enabled(self) -> ConfirmRunFinishResponse:
        # 1. 等待 transport 排空
        if self._transport is not None:
            if not self._transport.finish(timeout=None):
                return ConfirmRunFinishResponse(success=False, message="Transport is still running.")
            self._transport = None
        # 2. 停止心跳和上传跟踪器，释放资源
        if self._tracker is not None:
            self._tracker.set_state(CoreState.CORE_STATE_FINISHED)
        if self._heartbeat is not None:
            self._heartbeat.stop()
            self._heartbeat = None
        # 3. 上报最终的 finish_record 给后端，完成实验结束流程
        # 约定仅 online 模式暂存 finish_record，offline/local 模式在 deliver_run_finish 时就完成了全部流程，因此这里无需上报
        if self._mode != "online":
            return ConfirmRunFinishResponse(success=True, message="OK")
        if self._pending_online_finish_record is None:
            return ConfirmRunFinishResponse(
                success=False,
                message="Failed to confirm run finish: no pending finish record found.",
            )
        with safe.block(message="Failed to report run finish"):
            stop_experiment(
                self._ctx.username,
                self._ctx.project,
                self._ctx.experiment_id,
                state=self._pending_online_finish_record.state,
                finished_at=self._pending_online_finish_record.finished_at,
            )
            self._pending_online_finish_record = None
            return ConfirmRunFinishResponse(success=True, message="OK")
        # 虽然本地已经完成了全部流程，但由于网络等原因导致无法通知后端，因此返回失败状态，但是影响不大
        return ConfirmRunFinishResponse(success=False, message="Failed to finish run, but it has been saved locally.")
