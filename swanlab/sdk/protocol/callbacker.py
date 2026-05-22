"""
@author: cunyue
@file: callbacker.py
@time: 2026/4/12 01:12
@description: SwanLab SDK 回调器协议定义
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord

if TYPE_CHECKING:
    from swanlab.sdk.internal.settings import Settings


class Callback(ABC):
    """
    Base class for SwanLab callbacks, defining the lifecycle hooks of a SwanLab run.
    This class specifies the interface for custom callbacks that can be triggered
    at various stages of an experiment.
    """

    def on_run_initialized(self, run_dir: Path, path: str, settings: "Settings", *args, **kwargs) -> None:
        """
        Called immediately after ``swanlab.init`` has successfully executed.

        :param run_dir: The directory where the run is stored.
        :param path: The run path, formatted as ``/:username/:project/:run_id`` in online mode,
            or ``/:project/:run_id`` otherwise.
        """
        ...

    def on_scalar_flush(self, scalar_records: List[ScalarRecord], *args, **kwargs) -> None:
        """
        Called when new scalar records are flushed.

        :param scalar_records: A list of scalar records that have been flushed.
        """
        ...

    def on_media_flush(self, media_records: List[MediaRecord], *args, **kwargs) -> None:
        """
        Called when new media records are flushed.

        :param media_records: A list of media records that have been flushed.
        """
        ...

    def on_log_flush(self, log_records: List[LogRecord], *args, **kwargs) -> None:
        """
        Called when new log records are flushed.

        :param log_records: A list of log records that have been flushed.
        """
        ...

    def on_run_fork(self) -> None:
        """
        Called when the run is forked (e.g., in a multi-process setting).
        """
        ...

    def on_run_finished(self, state: str, error: Optional[str] = None, *args, **kwargs) -> None:
        """
        Called when the training or the SwanLab run finishes, either successfully or due to an error.

        :param state: The final state of the run (e.g., "finished", "crashed", "aborted").
        :param error: The error message or stack trace if the run exited abnormally, otherwise None.
        """
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """
        The globally unique identifier for this callback.

        This property is utilized by the callback operator for registration and deduplication.
        It must be implemented by all subclasses and return a unique string.

        :return: A unique string identifier.
        """
        ...
