"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13

@description: SwanLab Core Python 版本，封装SwanLab云端版核心业务，包括：
1. 提供http客户端，用于与SwanLab云端API进行交互。
2. 提供rpc封装函数，以rpc方式调用SwanLab云端API。
3. 提供上传线程，在另一个线程执行上传任务。
...

实现 CoreProtocol，当前为纯 Python 实现。
未来由 swanlab-core（Go 二进制）替代时，此模块整体被替换，
BackgroundConsumer 等调用方无需修改。
"""

from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import (
    FinishRecord,
    FinishResponse,
    StartRecord,
    StartResponse,
)
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core_python.api.experiment import create_or_resume_experiment, stop_experiment
from swanlab.sdk.internal.core_python.api.project import get_or_create_project, get_project
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.transport import Transport
from swanlab.sdk.internal.pkg import adapter, console, helper, safe
from swanlab.sdk.protocol import CoreProtocol
from swanlab.utils.experiment import generate_color, generate_name

__all__ = ["CorePython"]


class CorePython(CoreProtocol):
    """
    CoreProtocol 的 Python 实现。
    由 Run 在初始化时构造并注入给 BackgroundConsumer
    """

    START_RECORD_NUM = -1
    FINISH_RECORD_NUM = -2

    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        self._store: Optional[DataStoreWriter] = None
        self._transport: Optional[Transport] = None
        self._mode = ctx.config.settings.mode
        self._username: Optional[str] = None
        self._project: Optional[str] = None
        self._cuid: Optional[str] = None
        self._started: bool = False

    # ---------------------------------- start 方法 ----------------------------------

    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        if self._started:
            raise RuntimeError("Failed to start run: already started")
        resp = super().deliver_run_start(start_record)
        self._started = resp.success
        return resp

    def _start_store(self, resp: StartResponse):
        self._store = DataStoreWriter()
        self._store.open(str(self._ctx.run_file))
        record = Record(num=self.START_RECORD_NUM, start=resp.run)
        self._store.write(record.SerializeToString())

    def _start_without_cloud(self, start_record: StartRecord, message: str) -> StartResponse:
        resp = StartResponse(success=True, message=message, run=start_record)
        self._start_store(resp)
        return resp

    def _start_when_local(self, start_record: StartRecord) -> StartResponse:
        return self._start_without_cloud(start_record, "OK, but use local")

    def _start_when_offline(self, start_record: StartRecord) -> StartResponse:
        return self._start_without_cloud(start_record, "OK, but use offline")

    def _start_when_cloud(self, start_record: StartRecord) -> StartResponse:
        resp = self._report_run_start(start_record)
        self._start_store(resp)
        # Transport initialization is part of startup in cloud mode.
        # Fail fast on error instead of degrading silently.
        self._transport = Transport()
        return resp

    def _report_run_start(self, record: StartRecord) -> StartResponse:
        """
        运行开始
        :param record: 运行开始记录
        :return: 运行开始响应
        """
        # 1. 向后端同步运行开始事件
        # 获取当前项目，如果不存在则创建
        project_data = get_or_create_project(
            username=record.workspace,
            name=record.project,
            public=record.public,
        )
        username, project = project_data["username"], project_data["name"]
        # 获取当前项目详细信息
        project_info = get_project(username=username, name=project)
        # 获取当前实验
        history_experiment_count = project_info["_count"]["experiments"]
        name = record.name or generate_name(history_experiment_count)
        color = record.color or generate_color(history_experiment_count)
        resume = adapter.resume[record.resume]
        # 开启实验
        experiment = create_or_resume_experiment(
            username,
            project,
            name=name,
            resume=resume,
            run_id=record.id,
            color=color,
            description=record.description,
            job_type=record.job_type,
            group=record.group,
            tags=list(record.tags),
            created_at=record.started_at,
        )
        # 2. resume 时，向后端获取数据
        # TODO resume 时向后端获取数据或向本地获取数据

        # 3. 记录必要字段
        self._username = username
        self._project = project
        self._cuid = experiment["cuid"]
        # 4. 构建记录
        start_record = StartRecord()
        start_record.CopyFrom(record)
        start_record.name = name
        start_record.color = color
        start_record.resume = record.resume
        start_record.project = project
        start_record.workspace = username
        return StartResponse(success=True, message="OK", run=start_record)

    # ---------------------------------- publish 方法 ----------------------------------

    def publish(self, records: List[Record]) -> None:
        if not self._started:
            console.warning("CorePython is not started, skipping record publishing.")
            return
        super().publish(records)

    def _publish_store(self, records: List[Record]) -> None:
        assert self._store is not None, "store must be initialized before publishing"
        for record in records:
            self._store.write(record.SerializeToString())
            if helper.DEBUG:
                console.debug("Write record:", record.WhichOneof("record_type"))

    def _publish_when_local(self, records: List[Record]) -> None:
        self._publish_store(records)

    def _publish_when_offline(self, records: List[Record]) -> None:
        self._publish_store(records)

    def _publish_when_cloud(self, records: List[Record]) -> None:
        self._publish_store(records)
        assert self._transport is not None, "transport must be initialized before publishing"
        self._transport.put(records)

    # ---------------------------------- fork 方法 ----------------------------------

    def fork(self) -> "CorePython":
        raise RuntimeError("CorePython.fork() should not be called (designed for swanlab-core).")

    # ---------------------------------- finish 方法 ----------------------------------

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        if not self._started:
            raise RuntimeError("Failed to finish run: not started")
        resp = super().deliver_run_finish(finish_record)
        self._started = False
        return resp

    def _finish_store(self, record: Record):
        assert self._store is not None, "store must be initialized before shutdown"
        self._store.write(record.SerializeToString())
        self._store.close()
        self._store = None

    def _build_finish_record(self, finish_record: FinishRecord):
        record = Record(num=self.FINISH_RECORD_NUM)
        record.finish.CopyFrom(finish_record)
        return record

    def _finish_when_local(self, finish_record: FinishRecord) -> FinishResponse:
        record = self._build_finish_record(finish_record)
        self._finish_store(record)
        return FinishResponse(success=True, message="OK, but use local")

    def _finish_when_offline(self, finish_record: FinishRecord) -> FinishResponse:
        record = self._build_finish_record(finish_record)
        self._finish_store(record)
        return FinishResponse(success=True, message="OK, but use offline")

    def _finish_when_cloud(self, finish_record: FinishRecord) -> FinishResponse:
        record = self._build_finish_record(finish_record)
        self._finish_store(record)
        assert self._transport is not None, "transport must be initialized before finishing"
        self._transport.finish()
        self._transport = None
        resp = self._report_run_finish(finish_record)
        # 如果仅仅是与后端同步出现问题，则换一个让用户安心一些的提示信息
        if resp is None:
            return FinishResponse(success=False, message="Failed to finish run, but it has been saved locally.")
        return resp

    @safe.decorator(message="run finish error")
    def _report_run_finish(self, record: FinishRecord) -> FinishResponse:
        assert self._username is not None and self._project is not None and self._cuid is not None, (
            "Cannot finish cloud run: username, project, or cuid is missing."
        )
        # 向后端同步运行结束事件
        stop_experiment(
            self._username,
            self._project,
            self._cuid,
            state=record.state,
            finished_at=record.finished_at,
        )
        # 构建记录
        return FinishResponse(success=True, message="OK")
