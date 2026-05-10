"""
@author: cunyue
@file: core.py
@time: 2026/3/13
@description: SwanLab Core 接口协议，负责对产出的Record做持久化处理和后端交互
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import TYPE_CHECKING, List

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.env.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.run.v1.run_pb2 import (
    FinishRecord,
    FinishResponse,
    StartRecord,
    StartResponse,
)
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.pkg import safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext

__all__ = ["CoreEnum", "CoreProtocol"]


class CoreEnum(str, Enum):
    CORE_PYTHON = "CorePython"
    CORE = "Core"


class CoreProtocol(ABC):
    """
    SwanLab Core 协议
    协议层作为一个"垫片"，抹平前端与后端实现的差异，我们对CoreProtocol的错误处理理念基本与go一致，不raise Error，仅返回处理结果，上层根据返回结果做相应处理
    """

    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx
        self._mode = ctx.config.settings.mode

    # ---------------------------------- 实验开始 ----------------------------------

    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        """
        交付运行开始事件，同步事件，等待确认
        约定应该在此处完成Core相关组件的初始化，在这里完成不同模式的初始化函数分发
        """
        with safe.block(message="run start error"):
            if self._mode == "online":
                return self._start_when_online(start_record)
            elif self._mode == "local":
                return self._start_when_local(start_record)
            elif self._mode == "offline":
                return self._start_when_offline(start_record)
            return self._start_when_disabled(start_record)
        return StartResponse(success=False, message="Failed to start run")

    @staticmethod
    def _start_when_disabled(start_record: StartRecord) -> StartResponse:
        return StartResponse(success=True, message="I'm a teapot.", run=start_record)

    @abstractmethod
    def _start_when_local(self, start_record: StartRecord) -> StartResponse: ...

    @abstractmethod
    def _start_when_offline(self, start_record: StartRecord) -> StartResponse: ...

    @abstractmethod
    def _start_when_online(self, start_record: StartRecord) -> StartResponse: ...

    # ---------------------------------- 指标定义 ----------------------------------

    def upsert_columns(self, columns: List[ColumnRecord]) -> None:
        """上报一组 ColumnRecord，定义指标列"""
        with safe.block(message="upsert columns error"):
            if self._mode == "online":
                return self._upsert_columns_when_online(columns)
            elif self._mode == "local":
                return self._upsert_columns_when_local(columns)
            elif self._mode == "offline":
                return self._upsert_columns_when_offline(columns)
            return self._upsert_columns_when_disabled(columns)

    def _upsert_columns_when_disabled(self, columns: List[ColumnRecord]) -> None: ...

    @abstractmethod
    def _upsert_columns_when_local(self, columns: List[ColumnRecord]) -> None: ...

    @abstractmethod
    def _upsert_columns_when_offline(self, columns: List[ColumnRecord]) -> None: ...

    @abstractmethod
    def _upsert_columns_when_online(self, columns: List[ColumnRecord]) -> None: ...

    # ---------------------------------- 标量指标值 ----------------------------------

    def upsert_scalars(self, scalars: List[ScalarRecord]) -> None:
        """上报一组 ScalarRecord，记录标量指标值"""
        with safe.block(message="upsert data error"):
            if self._mode == "online":
                return self._upsert_scalars_when_online(scalars)
            elif self._mode == "local":
                return self._upsert_scalars_when_local(scalars)
            elif self._mode == "offline":
                return self._upsert_scalars_when_offline(scalars)
            return self._upsert_scalars_when_disabled(scalars)

    def _upsert_scalars_when_disabled(self, scalars: List[ScalarRecord]) -> None: ...

    @abstractmethod
    def _upsert_scalars_when_local(self, scalars: List[ScalarRecord]) -> None: ...

    @abstractmethod
    def _upsert_scalars_when_offline(self, scalars: List[ScalarRecord]) -> None: ...

    @abstractmethod
    def _upsert_scalars_when_online(self, scalars: List[ScalarRecord]) -> None: ...

    # ---------------------------------- 媒体指标值 ----------------------------------

    def upsert_media(self, media: List[MediaRecord]) -> None:
        """上报一组 MediaRecord，记录媒体指标值"""
        with safe.block(message="upsert data error"):
            if self._mode == "online":
                return self._upsert_media_when_online(media)
            elif self._mode == "local":
                return self._upsert_media_when_local(media)
            elif self._mode == "offline":
                return self._upsert_media_when_offline(media)
            return self._upsert_media_when_disabled(media)

    def _upsert_media_when_disabled(self, media: List[MediaRecord]) -> None: ...

    @abstractmethod
    def _upsert_media_when_local(self, media: List[MediaRecord]) -> None: ...

    @abstractmethod
    def _upsert_media_when_offline(self, media: List[MediaRecord]) -> None: ...

    @abstractmethod
    def _upsert_media_when_online(self, media: List[MediaRecord]) -> None: ...

    # ---------------------------------- 用户配置 ----------------------------------

    def upsert_configs(self, configs: List[ConfigRecord]) -> None:
        """上报一组 ConfigRecord，记录用户对 config 的修改"""
        with safe.block(message="upsert configs error"):
            if self._mode == "online":
                return self._upsert_configs_when_online(configs)
            elif self._mode == "local":
                return self._upsert_configs_when_local(configs)
            elif self._mode == "offline":
                return self._upsert_configs_when_offline(configs)
            return self._upsert_configs_when_disabled(configs)

    def _upsert_configs_when_disabled(self, configs: List[ConfigRecord]) -> None: ...

    @abstractmethod
    def _upsert_configs_when_local(self, configs: List[ConfigRecord]) -> None: ...

    @abstractmethod
    def _upsert_configs_when_offline(self, configs: List[ConfigRecord]) -> None: ...

    @abstractmethod
    def _upsert_configs_when_online(self, configs: List[ConfigRecord]) -> None: ...

    # ---------------------------------- 终端输出 ----------------------------------

    def upsert_logs(self, logs: List[LogRecord]) -> None:
        """上报一组 LogRecord，记录用户终端输出"""
        with safe.block(message="upsert logs error"):
            if self._mode == "online":
                return self._upsert_logs_when_online(logs)
            elif self._mode == "local":
                return self._upsert_logs_when_local(logs)
            elif self._mode == "offline":
                return self._upsert_logs_when_offline(logs)
            return self._upsert_logs_when_disabled(logs)

    def _upsert_logs_when_disabled(self, logs: List[LogRecord]) -> None: ...

    @abstractmethod
    def _upsert_logs_when_local(self, logs: List[LogRecord]) -> None: ...

    @abstractmethod
    def _upsert_logs_when_offline(self, logs: List[LogRecord]) -> None: ...

    @abstractmethod
    def _upsert_logs_when_online(self, logs: List[LogRecord]) -> None: ...

    # ---------------------------------- 依赖信息 ----------------------------------

    def upsert_requirements(self, requirements: List[RequirementsRecord]) -> None:
        """上报一组 RequirementsRecord，记录对 requirements 的修改"""
        with safe.block(message="upsert requirements error"):
            if self._mode == "online":
                return self._upsert_requirements_when_online(requirements)
            elif self._mode == "local":
                return self._upsert_requirements_when_local(requirements)
            elif self._mode == "offline":
                return self._upsert_requirements_when_offline(requirements)
            return self._upsert_requirements_when_disabled(requirements)

    def _upsert_requirements_when_disabled(self, requirements: List[RequirementsRecord]) -> None: ...

    @abstractmethod
    def _upsert_requirements_when_local(self, requirements: List[RequirementsRecord]) -> None: ...

    @abstractmethod
    def _upsert_requirements_when_offline(self, requirements: List[RequirementsRecord]) -> None: ...

    @abstractmethod
    def _upsert_requirements_when_online(self, requirements: List[RequirementsRecord]) -> None: ...

    # ---------------------------------- conda信息 ----------------------------------

    def upsert_conda(self, conda: List[CondaRecord]) -> None:
        """上报一组 CondaRecord，记录对 conda 的修改"""
        with safe.block(message="upsert conda error"):
            if self._mode == "online":
                return self._upsert_conda_when_online(conda)
            elif self._mode == "local":
                return self._upsert_conda_when_local(conda)
            elif self._mode == "offline":
                return self._upsert_conda_when_offline(conda)
            return self._upsert_conda_when_disabled(conda)

    def _upsert_conda_when_disabled(self, conda: List[CondaRecord]) -> None: ...

    @abstractmethod
    def _upsert_conda_when_local(self, conda: List[CondaRecord]) -> None: ...

    @abstractmethod
    def _upsert_conda_when_offline(self, conda: List[CondaRecord]) -> None: ...

    @abstractmethod
    def _upsert_conda_when_online(self, conda: List[CondaRecord]) -> None: ...

    # ---------------------------------- 系统元信息 ----------------------------------

    def upsert_metadata(self, metadata: List[MetadataRecord]) -> None:
        """上报一组 MetadataRecord，记录对 metadata 的修改"""
        with safe.block(message="upsert metadata error"):
            if self._mode == "online":
                return self._upsert_metadata_when_online(metadata)
            elif self._mode == "local":
                return self._upsert_metadata_when_local(metadata)
            elif self._mode == "offline":
                return self._upsert_metadata_when_offline(metadata)
            return self._upsert_metadata_when_disabled(metadata)

    def _upsert_metadata_when_disabled(self, metadata: List[MetadataRecord]) -> None: ...

    @abstractmethod
    def _upsert_metadata_when_local(self, metadata: List[MetadataRecord]) -> None: ...

    @abstractmethod
    def _upsert_metadata_when_offline(self, metadata: List[MetadataRecord]) -> None: ...

    @abstractmethod
    def _upsert_metadata_when_online(self, metadata: List[MetadataRecord]) -> None: ...

    # ---------------------------------- 文件保存 ----------------------------------

    def upsert_saves(self, saves: List[SaveRecord]) -> None:
        """上报一组 SaveRecord，记录用户通过 swanlab.save() 保存的文件"""
        with safe.block(message="upsert saves error"):
            if self._mode == "online":
                return self._upsert_saves_when_online(saves)
            elif self._mode == "local":
                return self._upsert_saves_when_local(saves)
            elif self._mode == "offline":
                return self._upsert_saves_when_offline(saves)
            return self._upsert_saves_when_disabled(saves)

    def _upsert_saves_when_disabled(self, saves: List[SaveRecord]) -> None: ...

    @abstractmethod
    def _upsert_saves_when_local(self, saves: List[SaveRecord]) -> None: ...

    @abstractmethod
    def _upsert_saves_when_offline(self, saves: List[SaveRecord]) -> None: ...

    @abstractmethod
    def _upsert_saves_when_online(self, saves: List[SaveRecord]) -> None: ...

    # ---------------------------------- 运行结束 ----------------------------------

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        """
        交付运行结束事件，同步事件，等待确认
        约定应该在此处完成Core相关组件的清理，在这里完成不同模式的结束函数分发
        """
        with safe.block(message="run finish error"):
            if self._mode == "online":
                return self._finish_when_online(finish_record)
            elif self._mode == "local":
                return self._finish_when_local(finish_record)
            elif self._mode == "offline":
                return self._finish_when_offline(finish_record)
            return self._finish_when_disabled()
        return FinishResponse(success=False, message="Failed to finish run")

    @staticmethod
    def _finish_when_disabled() -> FinishResponse:
        return FinishResponse(success=True, message="I'm a teapot.")

    @abstractmethod
    def _finish_when_local(self, finish_record: FinishRecord) -> FinishResponse: ...

    @abstractmethod
    def _finish_when_offline(self, finish_record: FinishRecord) -> FinishResponse: ...

    @abstractmethod
    def _finish_when_online(self, finish_record: FinishRecord) -> FinishResponse: ...

    # ---------------------------------- 进程fork ----------------------------------

    @abstractmethod
    def fork(self) -> "CoreProtocol":
        """
        创建一个新的 CoreProtocol，此函数用于处理多进程fork的情况
        """
        ...
