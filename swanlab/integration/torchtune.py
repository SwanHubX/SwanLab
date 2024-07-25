from torchtune.utils._distributed import get_world_size_and_rank
from torchtune.utils.metric_logging import MetricLoggerInterface
from typing import Mapping, Optional, Union

from numpy import ndarray
from omegaconf import DictConfig, OmegaConf
from torch import Tensor

Scalar = Union[Tensor, ndarray, int, float]


class SwanLabLogger(MetricLoggerInterface):
    """
    Example:
        In torchtune config:
        >>> # Logging
        >>> metric_logger:
        >>>     _component_: swanlab.integration.huggingface.SwanLabLogger
        >>>     project: "gemma-lora-finetune"
        >>>     experiment_name: "gemma-2b"
        >>>     log_dir: ${output_dir}

    Raises:
        ImportError: If ``swanlab`` package is not installed.

    Note:
        This logger requires the swanlab package to be installed.
        You can install it with `pip install swanlab`.
        In order to use the logger, you need to login to your SwanLab account.
        You can do this by running `swanlab login` in your terminal.
    """

    def __init__(
        self,
        project: str = "torchtune",
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        mode: Optional[str] = None,
        log_dir: Optional[str] = None,
        **kwargs,
    ):
        try:
            import swanlab
        except ImportError as e:
            raise ImportError(
                "``swanlab`` package not found. Please install swanlab using `pip install swanlab` to use SwanLabLogger."
            ) from e
        self._swanlab = swanlab

        # Use dir if specified, otherwise use log_dir.
        self.log_dir = kwargs.pop("dir", log_dir)

        _, self.rank = get_world_size_and_rank()

        if self._swanlab.get_run() is None and self.rank == 0:
            run = self._swanlab.init(
                project=project,
                workspace=workspace,
                name=experiment_name,
                description=description,
                mode=mode,
                logdir=self.log_dir,
                **kwargs,
            )

        self.config_allow_val_change = kwargs.get("allow_val_change", False)

    def log_config(self, config: DictConfig) -> None:
        if self._swanlab.get_run():
            resolved = OmegaConf.to_container(config, resolve=True)
            self._swanlab.config.update(resolved, allow_val_change=self.config_allow_val_change)

    def log(self, name: str, data: Scalar, step: int) -> None:
        if self._swanlab.get_run():
            self._swanlab.log({name: data}, step=step)

    def log_dict(self, payload: Mapping[str, Scalar], step: int) -> None:
        if self._swanlab.get_run():
            self._swanlab.log({**payload}, step=step)

    def __del__(self) -> None:
        if self._swanlab.get_run():
            self._swanlab.finish()

    def close(self) -> None:
        if self._swanlab.get_run():
            self._swanlab.finish()
