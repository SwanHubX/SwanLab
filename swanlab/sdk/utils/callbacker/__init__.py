"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 17:30
@description: SwanLab 运行时回调函数
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

__all__ = ["SwanLabCallback"]


class SwanLabCallback(ABC):
    """
    Base class for SwanLab callbacks, defining the lifecycle hooks of a SwanLab run.
    This class specifies the interface for custom callbacks that can be triggered
    at various stages of an experiment.
    """

    def on_run_init(self, logdir: str, path: str) -> None:
        """
        Called immediately after `swanlab.init` has successfully executed.

        :param logdir: The directory path where the user's local logs are stored.
        :param path: The routing path of the experiment, formatted as `/:workspace/:project/:run_id`.
        """
        pass

    def on_metadata_update(self, *args, **kwargs) -> None:
        """
        Called when the experiment's metadata is updated.
        This includes updates to:
        - Experiment configuration (RunConfig)
        - Hardware/Host information
        - Environment information
        """
        pass

    def on_log(self, data: Dict[str, Any], step: Optional[int], *args, **kwargs) -> None:
        """
        Called every time `swanlab.log` is executed to record metrics or media.

        :param data: The payload dictionary containing the logged key-value pairs.
        :param step: The global step at which the data was logged. Can be None if not explicitly tracked.
        """
        pass

    def on_column_created(self, column: Any, *args, **kwargs) -> None:
        """
        Called when a new column is added to the experiment tracking schema.

        :param column: The column object containing the newly created schema information.
        """
        pass

    def on_metric_created(self, metric: Any, *args, **kwargs) -> None:
        """
        Called when a new metric tracker is initialized.

        :param metric: The metric object representing the newly created metric.
        """
        pass

    def on_run_end(self, state: str, error: Optional[str] = None) -> None:
        """
        Called when the training or the SwanLab run finishes, either successfully or due to an error.

        :param state: The final state of the run (e.g., "finished", "crashed", "aborted").
        :param error: The error message or stack trace if the run exited abnormally, otherwise None.
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """
        The globally unique identifier for this callback.

        This property is utilized by the callback operator for registration and deduplication.
        It must be implemented by all subclasses and return a unique string.

        :return: A unique string identifier.
        """
        pass
