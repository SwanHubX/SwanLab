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

from typing import Callable, List, Optional

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
from swanlab.sdk.internal.pkg.safe import block as safe_block
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

    def __init__(self, ctx: RunContext, upload_callback: Optional[Callable[[int], None]] = None):
        super().__init__(ctx)
        self._store: Optional[DataStoreWriter] = None
        self._transport: Optional[Transport] = None
        self._mode = ctx.config.settings.mode
        self._upload_callback = upload_callback
        self._username: Optional[str] = None
        self._project: Optional[str] = None
        self._cuid: Optional[str] = None

    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        if self._store is not None or self._transport is not None:
            raise RuntimeError("CorePython has already been started.")
        # 1. 向后端同步运行开始事件
        resp = self._run_start(start_record)
        if resp is None:
            return StartResponse(success=False, message="Failed to start run.")
        # 2. 启动组件
        if self._mode != "disabled":
            self._store = DataStoreWriter()
            self._store.open(str(self._ctx.run_file))
            record = Record(num=self.START_RECORD_NUM, start=resp.run)
            self._store.write(record.SerializeToString())
        if self._mode == "cloud":
            self._transport = Transport(upload_callback=self._upload_callback)
        return resp

    def publish(self, records: List[Record]) -> None:
        if self._store is None and self._transport is None:
            console.warning("CorePython is not started, skipping record handling.")
            return
        with safe_block(message="CorePython publish error"):
            if self._store is not None:
                for record in records:
                    self._store.write(record.SerializeToString())
                    if helper.DEBUG:
                        console.debug("Write record:", record.WhichOneof("record_type"))
            if self._transport is not None:
                self._transport.put(records)

    def fork(self) -> "CorePython":
        raise NotImplementedError(
            "CorePython.fork() is not implemented. Please waiting for go version, while you should not reach here?"
        )

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        # 1. 构建记录
        record = Record(
            num=self.FINISH_RECORD_NUM,
            finish=FinishRecord(
                state=finish_record.state, error=finish_record.error, finished_at=finish_record.finished_at
            ),
        )
        # 2. 停止组件
        if self._transport is not None:
            self._transport.finish()
            self._transport = None
        use_store = self._store is not None
        if self._store is not None:
            self._store.write(record.SerializeToString())
            self._store.close()
            self._store = None
        # 3. 向后端同步运行结束事件
        result = self._run_finish(finish_record)
        if result is None:
            # 如果与后端同步失败，根据不同情况，返回不同的响应
            if use_store:
                return FinishResponse(
                    success=False,
                    message="Failed to sync the run finish event to the server, but it has been saved locally.",
                )

            else:
                return FinishResponse(
                    success=False,
                    message="Failed to sync the run finish event to the server.",
                )
        return FinishResponse(success=True, message="Done")

    @safe.decorator(message="run start error")
    def _run_start(self, record: StartRecord) -> StartResponse:
        """
        运行开始
        :param record: 运行开始记录
        :return: 运行开始响应
        """
        # 0. 如果不是 cloud 模式则直接返回
        if self._mode != "cloud":
            return StartResponse(success=True, message="OK, but no cloud", run=record)
        # 1. 向后端同步运行开始事件
        # 获取当前项目，如果不存在则创建
        project = get_or_create_project(
            username=record.workspace,
            name=record.project,
            public=record.public,
        )
        username, project = project["username"], project["name"]
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

    @safe.decorator(message="run finish error")
    def _run_finish(self, record: FinishRecord) -> FinishResponse:
        """
        运行结束
        :param record: 运行结束请求
        :return: 运行结束响应
        """
        # 如果不是 cloud 模式则直接返回
        if self._mode != "cloud":
            return FinishResponse(success=True, message="OK, but no cloud")

        assert self._username is not None and self._project is not None and self._cuid is not None, (
            "Required fields are not set."
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
