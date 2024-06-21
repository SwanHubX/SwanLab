"""
swanlab.integration.torchtune.SwanLabLogger
"""

from typing import Mapping, Optional
from omegaconf import DictConfig, OmegaConf
from torchtune.utils.metric_logging import MetricLoggerInterface, Scalar, get_world_size_and_rank
import swanlab


class SwanLabLogger(MetricLoggerInterface):
    def __init__(
        self,
        project: str = "torchtune",
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        mode: Optional[str] = None,
        logdir: Optional[str] = None,
        **kwargs,
    ):

        self._swanlab = swanlab

        _, self.rank = get_world_size_and_rank()

        if self._swanlab.get_run() is None and self.rank == 0:
            self._swanlab.init(
                project=project,
                workspace=workspace,
                experiment_name=experiment_name,
                description=description,
                logdir=logdir,
                mode=mode,
                **kwargs,
            )

    def log_config(self, config: DictConfig) -> None:
        if self._swanlab.get_run():
            resolved = OmegaConf.to_container(config, resolve=True)
            self._swanlab.config.update(resolved)

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
