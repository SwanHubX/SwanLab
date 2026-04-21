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
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.run.v1.run_pb2 import (
    FinishRecord,
    FinishResponse,
    StartRecord,
    StartResponse,
)
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.proto.swanlab.system.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
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
            if self._mode == "cloud":
                return self._start_when_cloud(start_record)
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
    def _start_when_cloud(self, start_record: StartRecord) -> StartResponse: ...

    # ---------------------------------- 指标定义 ----------------------------------

    def upsert_columns(self, columns: List[ColumnRecord]) -> None:
        """上报一组 ColumnRecord，定义指标列"""
        with safe.block(message="upsert columns error"):
            if self._mode == "cloud":
                return self._upsert_columns_when_cloud(columns)
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
    def _upsert_columns_when_cloud(self, columns: List[ColumnRecord]) -> None: ...

    # ---------------------------------- 指标值 ----------------------------------

    def upsert_data(self, data: List[DataRecord]) -> None:
        """上报一组 DataRecord，记录指标值"""
        with safe.block(message="upsert data error"):
            if self._mode == "cloud":
                return self._upsert_data_when_cloud(data)
            elif self._mode == "local":
                return self._upsert_data_when_local(data)
            elif self._mode == "offline":
                return self._upsert_data_when_offline(data)
            return self._upsert_data_when_disabled(data)

    def _upsert_data_when_disabled(self, data: List[DataRecord]) -> None: ...

    @abstractmethod
    def _upsert_data_when_local(self, data: List[DataRecord]) -> None: ...

    @abstractmethod
    def _upsert_data_when_offline(self, data: List[DataRecord]) -> None: ...

    @abstractmethod
    def _upsert_data_when_cloud(self, data: List[DataRecord]) -> None: ...

    # ---------------------------------- 用户配置 ----------------------------------

    def upsert_configs(self, configs: List[ConfigRecord]) -> None:
        """上报一组 ConfigRecord，记录用户对 config 的修改"""
        with safe.block(message="upsert configs error"):
            if self._mode == "cloud":
                return self._upsert_configs_when_cloud(configs)
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
    def _upsert_configs_when_cloud(self, configs: List[ConfigRecord]) -> None: ...

    # ---------------------------------- 终端输出 ----------------------------------

    def upsert_consoles(self, consoles: List[ConsoleRecord]) -> None:
        """上报一组 ConsoleRecord，记录用户终端输出"""
        with safe.block(message="upsert consoles error"):
            if self._mode == "cloud":
                return self._upsert_consoles_when_cloud(consoles)
            elif self._mode == "local":
                return self._upsert_consoles_when_local(consoles)
            elif self._mode == "offline":
                return self._upsert_consoles_when_offline(consoles)
            return self._upsert_consoles_when_disabled(consoles)

    def _upsert_consoles_when_disabled(self, consoles: List[ConsoleRecord]) -> None: ...

    @abstractmethod
    def _upsert_consoles_when_local(self, consoles: List[ConsoleRecord]) -> None: ...

    @abstractmethod
    def _upsert_consoles_when_offline(self, consoles: List[ConsoleRecord]) -> None: ...

    @abstractmethod
    def _upsert_consoles_when_cloud(self, consoles: List[ConsoleRecord]) -> None: ...

    # ---------------------------------- 依赖信息 ----------------------------------

    def upsert_requirements(self, requirements: List[RequirementsRecord]) -> None:
        """上报一组 RequirementsRecord，记录对 requirements 的修改"""
        with safe.block(message="upsert requirements error"):
            if self._mode == "cloud":
                return self._upsert_requirements_when_cloud(requirements)
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
    def _upsert_requirements_when_cloud(self, requirements: List[RequirementsRecord]) -> None: ...

    # ---------------------------------- conda信息 ----------------------------------

    def upsert_conda(self, conda: List[CondaRecord]) -> None:
        """上报一组 CondaRecord，记录对 conda 的修改"""
        with safe.block(message="upsert conda error"):
            if self._mode == "cloud":
                return self._upsert_conda_when_cloud(conda)
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
    def _upsert_conda_when_cloud(self, conda: List[CondaRecord]) -> None: ...

    # ---------------------------------- 系统元信息 ----------------------------------

    def upsert_metadata(self, metadata: List[MetadataRecord]) -> None:
        """上报一组 MetadataRecord，记录对 metadata 的修改"""
        with safe.block(message="upsert metadata error"):
            if self._mode == "cloud":
                return self._upsert_metadata_when_cloud(metadata)
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
    def _upsert_metadata_when_cloud(self, metadata: List[MetadataRecord]) -> None: ...

    # ---------------------------------- 运行结束 ----------------------------------

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        """
        交付运行结束事件，同步事件，等待确认
        约定应该在此处完成Core相关组件的清理，在这里完成不同模式的结束函数分发
        """
        with safe.block(message="run finish error"):
            if self._mode == "cloud":
                return self._finish_when_cloud(finish_record)
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
    def _finish_when_cloud(self, finish_record: FinishRecord) -> FinishResponse: ...

    # ---------------------------------- 进程fork ----------------------------------

    @abstractmethod
    def fork(self) -> "CoreProtocol":
        """
        创建一个新的 CoreProtocol，此函数用于处理多进程fork的情况
        """
        ...
