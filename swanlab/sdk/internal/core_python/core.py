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
from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord, SaveType
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.api.experiment import (
    get_experiment_summary,
    stop_experiment,
)
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.heartbeat import Heartbeat
from swanlab.sdk.internal.core_python.metrics import RunMetrics
from swanlab.sdk.internal.core_python.pkg import builder, counter
from swanlab.sdk.internal.core_python.store import DataStoreWriter, DataStoreWriterProtocol, NullDataStoreWriter
from swanlab.sdk.internal.core_python.transport import Transport
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.utils import generate_run_online_path, prepare_experiment_start
from swanlab.sdk.internal.core_python.watcher import (
    FileWatcher,
    FileWatcherProtocol,
    NullFileWatcher,
    create_save_links,
)
from swanlab.sdk.internal.pkg import adapter, console, safe
from swanlab.sdk.internal.pkg.client import Client
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
        self._run_ctx: Optional[CoreContext] = None
        self._store: Optional[DataStoreWriterProtocol] = None
        self._transport: Optional[Transport] = None
        # 标记core是否激活，未激活时拒绝接受上报数据，表达一个完整的生命周期状态：
        # initialized but not started -> active -> finished
        self._active = False
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
        self._watcher: Optional[FileWatcherProtocol] = None
        self._pending_end_saves: List[Record] = []
        # 按实例身份释放本 Core 创建的 client。
        self._client_instance: Optional[Client] = None

    @property
    def _ctx(self) -> CoreContext:
        assert self._run_ctx, "run context not set"
        return self._run_ctx

    @_ctx.setter
    def _ctx(self, ctx: CoreContext):
        self._run_ctx = ctx

    # ---------------------------------- 实验开始 ----------------------------------

    def deliver_run_start(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        if self._active:
            return DeliverRunStartResponse(success=False, message="Run has already been active.")
        resp = super().deliver_run_start(start_request)
        self._active = resp.success
        return resp

    def _start_store(self, resp: DeliverRunStartResponse):
        # skip_store 模式下使用空写入器与空监听器，跳过本地存储和文件监听
        if self._ctx.config.skip_store:
            self._store = NullDataStoreWriter()
            self._watcher = NullFileWatcher()
        else:
            self._store = DataStoreWriter()
            self._watcher = FileWatcher(on_change=self._on_file_changed)
        self._store.open(str(self._ctx.run_file))
        self._store_records([builder.build_start_record(resp.run)])

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
        # 1. 进程内独占互斥：联网 Core（Run/Sync）至多一个，client 单例即原子占用位。
        #    正常流程下此处不应存在残留单例（init 不再预建、login 不驻留、sync 用后即毁），
        #    该检查用于在创建任何资源前给出可读的拒绝信息，而非撞上 client.new 的 already exists。
        if client.exists():
            return DeliverRunStartResponse(
                success=False,
                message=(
                    "A networked SwanLab core is already active in this process. Finish it before starting a new one."
                ),
            )
        # 2. 凭证合法性校验：None / 空字符串 / 仅空白均视为无效，直接阻断启动
        api_key = self._ctx.config.api_key
        if api_key is None or not api_key.strip():
            return DeliverRunStartResponse(
                success=False,
                message=(
                    "Online mode requires a valid API key. "
                    "Login via `swanlab.login()` or configure SWANLAB_API_KEY / .netrc credentials."
                ),
            )
        # 3. 自建 client：仅由 Core 基于 Proto 下沉的凭据创建，创建成功后记录所有权
        try:
            owned_instance = client.new(api_key=api_key, base_url=self._ctx.config.api_host)
            self._client_instance = owned_instance
            resp = self._report_run_start(start_request.start_record)
            self._start_store(resp)
            self._tracker = UploadTracker()
            self._tracker.set_state(CoreState.CORE_STATE_RUNNING)
            self._transport = Transport(ctx=self._ctx, tracker=self._tracker)
            self._heartbeat = Heartbeat(self._ctx.experiment_id)
            self._heartbeat.start()
            return resp
        except BaseException:
            # 启动半途失败必须严格回滚：先停止已启动的消费者，再集中释放单例，严禁残留活动单例
            self._rollback_run_start()
            raise

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
            project_id=run_info.project_data["cuid"],
            project_version=run_info.project_data.get("version", None),
            experiment_id=run_info.experiment["cuid"],
        )
        # 3. resume 时，向后端获取数据
        summary: Optional[ResumeExperimentSummaryType] = None
        if not run_info.new_experiment:
            summary = get_experiment_summary(
                self._ctx.project_id,
                self._ctx.experiment_id,
                created_at=run_info.experiment["createdAt"],
            )
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

    # ---------------------------------- client 生命周期 ----------------------------------

    def _rollback_run_start(self) -> None:
        """
        启动失败回滚：停止已启动的在线消费者、关闭本地存储，最后集中释放 client 单例。
        严格确认消费者退出后才允许释放：未确认退出时保留所有权，不允许新联网 Core 启动。
        """
        consumers_stopped = self._shutdown_online_consumers(timeout=Transport.FINISH_JOIN_TIMEOUT)
        if self._watcher is not None:
            with safe.block(message="Failed to stop file watcher during rollback", level="debug"):
                self._watcher.stop()
            self._watcher = None
        if self._store is not None:
            with safe.block(message="Failed to close local store during rollback", level="debug"):
                self._store.close()
            self._store = None
        # 仅当消费者确认退出时才允许释放 client；未确认时保留占用和所有权
        if consumers_stopped:
            self._teardown_client()
        else:
            console.warning(
                "Run rollback could not confirm all consumers stopped. "
                "Client ownership retained to prevent concurrent Core startup with active workers."
            )

    def _shutdown_online_consumers(self, timeout: Optional[float]) -> bool:
        """Stop consumers, retaining references and ownership unless every stop succeeds."""
        if self._transport is not None:
            with safe.block(message="Failed to stop transport during shutdown"):
                if self._transport.finish(timeout=timeout):
                    self._transport = None
        if self._transport is not None:
            return False
        if self._heartbeat is not None:
            with safe.block(message="Failed to stop heartbeat during shutdown"):
                self._heartbeat.stop()
                self._heartbeat = None
        stopped = self._transport is None and self._heartbeat is None
        if stopped and self._tracker is not None:
            self._tracker.set_state(CoreState.CORE_STATE_FINISHED)
        return stopped

    def _teardown_client(self) -> bool:
        if self._client_instance is None:
            return True
        with safe.block(message="Failed to release SwanLab client"):
            client.reset(client=self._client_instance)
            self._client_instance = None
        return self._client_instance is None

    # ---------------------------------- 数据上报 ----------------------------------

    def _store_records(self, records: List[Record]) -> None:
        """将一组 Record 写入本地存储。"""
        assert self._store is not None, "store must be initialized before upsert"
        if self._ctx.config.skip_store:
            # skip_store：跳过本地序列化与写入
            return
        for record in records:
            self._store.write(record.SerializeToString())

    def _transport_put(self, records: List[Record]) -> None:
        """将一组 Record 推送到上传队列"""
        assert self._transport is not None, "transport must be initialized before upsert"
        self._transport.put(records)

    # ---- upsert_columns ----

    def upsert_columns(self, columns: List[ColumnRecord]) -> None:
        if not self._active:
            console.warning("Core is not active, refusing to upsert columns")
            return
        return super().upsert_columns(columns)

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

    def upsert_scalars(self, scalars: List[ScalarRecord]) -> None:
        if not self._active:
            console.warning("Core is not active, refusing to upsert scalars")
            return
        return super().upsert_scalars(scalars)

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

    def upsert_media(self, media: List[MediaRecord]) -> None:
        if not self._active:
            console.warning("Core is not active, refusing to upsert media")
            return
        return super().upsert_media(media)

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

    def upsert_logs(self, logs: List[LogRecord]) -> None:
        if not self._active:
            console.warning("Core is not active, refusing to upsert logs")
            return
        return super().upsert_logs(logs)

    def _upsert_logs_when_local(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)

    def _upsert_logs_when_offline(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)

    def _upsert_logs_when_online(self, logs: List[LogRecord]) -> None:
        records = [builder.build_log_record(self._counter, self._epoch, c) for c in logs]
        self._store_records(records)
        if records:
            action = "Skipped storing" if self._ctx.config.skip_store else "Stored"
            console.debug(
                f"{action} log records locally: count={len(records)}, nums={records[0].num}..{records[-1].num}",
                write_to_tty=False,
            )
        self._transport_put(records)
        if records:
            console.debug(
                f"Queued log records for upload: count={len(records)}, nums={records[0].num}..{records[-1].num}",
                write_to_tty=False,
            )

    # ---- upsert_saves ----

    def upsert_saves(self, saves: List[SaveRecord]) -> None:
        if not self._active:
            console.warning("Core is not active, refusing to upsert saves")
            return
        if self._mode == "disabled":
            return self._upsert_saves_when_disabled(saves)
        # 内部保存（如config、metadata等）和用户保存（文件保存）共用 SaveRecord 结构
        # 目前它们在产品设计上暂未统一，换句话说上传config、metadata文件的时候并不会保存对应文件
        # 因此需要区分两者，内部保存仅上传，用户保存走save逻辑
        # 1. 区分内部保存和用户保存
        custom_saves: List[SaveRecord] = []
        internal_saves: List[SaveRecord] = []
        for save in saves:
            if save.type == SaveType.SAVE_TYPE_CUSTOM:
                custom_saves.append(save)
            else:
                internal_saves.append(save)
        # 2. 执行内部保存的上传逻辑：本地持久化存储、上传（如果在线）等，但不注册文件监视器
        if internal_saves:
            records = [builder.build_save_record(self._counter, s, s.type) for s in internal_saves]
            self._store_records(records)
            if self._mode == "online":
                self._transport_put(records)
        # 3. 执行用户保存的逻辑
        if custom_saves:
            if self._mode == "local":
                self._upsert_saves_when_local(custom_saves)
            elif self._mode == "offline":
                self._upsert_saves_when_offline(custom_saves)
            elif self._mode == "online":
                self._upsert_saves_when_online(custom_saves)

    def _handle_custom_save(self, saves: List[SaveRecord]) -> List[Record]:
        # skip_store：不创建本地镜像软链接（不触碰 files_dir、不填 target_path），无文件监听
        if not self._ctx.config.skip_store:
            assert self._watcher is not None, "watcher must be initialized before upsert"
            linked = create_save_links(saves, self._ctx.files_dir)
            if linked > 0:
                console.info(
                    f"Symlinked {linked} files into the SwanLab run directory; call swanlab.save again to sync new files."
                )
            self._watcher.register_live_watches(saves, self._ctx.files_dir)
        records = [builder.build_save_record(self._counter, s) for s in saves]
        self._store_records(records)
        return records

    def _upsert_saves_when_local(self, saves: List[SaveRecord]) -> None:
        self._handle_custom_save(saves)

    def _upsert_saves_when_offline(self, saves: List[SaveRecord]) -> None:
        self._handle_custom_save(saves)

    def _upsert_saves_when_online(self, saves: List[SaveRecord]) -> None:
        records = self._handle_custom_save(saves)
        # "end" policy saves: store locally but defer cloud upload until finish
        transport_records: List[Record] = []
        for save, record in zip(saves, records):
            if save.policy == SavePolicy.SAVE_POLICY_END:
                self._pending_end_saves.append(record)
            else:
                transport_records.append(record)
        if transport_records:
            self._transport_put(transport_records)

    def _on_file_changed(self, save_record: SaveRecord) -> None:
        """FileWatcher 回调：持久化 + 可选上传。"""
        if not self._active or self._store is None:
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
        if not self._active:
            return DeliverRunFinishResponse(success=False, message="Run is not active, refusing to finish.")
        resp = super().deliver_run_finish(finish_request)
        self._active = False
        if not resp.success:
            self._rollback_run_start()
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
        self._store_records([record])
        self._store.close()
        self._store = None
        if self._watcher is not None:
            self._watcher.stop()
            self._watcher = None
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

    def _get_operation_when_online(self) -> GetOperationStatsResponse:
        assert self._tracker is not None, "tracker must be initialized before get_operation_stats"
        stats = self._tracker.snapshot()
        return GetOperationStatsResponse(
            success=True,
            message="OK",
            stats=stats,
        )

    # ---------------------------------- 确认完成 ----------------------------------

    def _confirm_finish_when_enabled(self) -> ConfirmRunFinishResponse:
        if self._active:
            return ConfirmRunFinishResponse(success=False, message="Run is still active, cannot confirm finish.")

        response = ConfirmRunFinishResponse(success=False, message="Transport or heartbeat is still running.")
        consumers_stopped = False
        try:
            consumers_stopped = self._shutdown_online_consumers(timeout=None)
            if not consumers_stopped:
                return response
            # 3. 上报最终的 finish_record 给后端，完成实验结束流程
            # 约定仅 online 模式暂存 finish_record，offline/local 模式在 deliver_run_finish 时就完成了全部流程，因此这里无需上报
            if self._mode != "online":
                response.success = True
                response.message = "OK"
                return response
            if self._pending_online_finish_record is None:
                response.message = "Failed to confirm run finish: no pending finish record found."
                return response
            with safe.block(message="Failed to report run finish"):
                stop_experiment(
                    self._ctx.username,
                    self._ctx.project,
                    self._ctx.experiment_id,
                    state=self._pending_online_finish_record.state,
                    finished_at=self._pending_online_finish_record.finished_at,
                )
                self._pending_online_finish_record = None
                response.success = True
                response.message = "OK"
                return response
            # 虽然本地已经完成了全部流程，但由于网络等原因导致无法通知后端，因此返回失败状态，但是影响不大。
            # skip_store 下没有本地副本，提示语不能暗示数据仍可从本地恢复
            if self._ctx.config.skip_store:
                message = (
                    "Failed to finish run, and no local copy was kept (skip_store); "
                    "the run state on the cloud is unconfirmed."
                )
            else:
                message = "Failed to finish run, but it has been saved locally."
            response.message = message
            return response
        except BaseException:
            consumers_stopped = self._shutdown_online_consumers(timeout=Transport.FINISH_JOIN_TIMEOUT)
            if not consumers_stopped:
                console.warning(
                    "Run finish interrupted: could not confirm all consumers stopped. Client ownership retained."
                )
            raise
        finally:
            if consumers_stopped and not self._teardown_client():
                response.success = False
                response.message = "Failed to release SwanLab client; the core still owns its client."
