from pathlib import Path
from typing import Mapping, Optional
from omegaconf import DictConfig, OmegaConf
from torchtune.utils.metric_logging import MetricLoggerInterface, Scalar, get_world_size_and_rank
import swanlab
import os


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

        # Use dir if specified, otherwise use log_dir.
        self.log_dir = kwargs.pop("dir", logdir)

        _, self.rank = get_world_size_and_rank()

        if self._swanlab.get_run() is None and self.rank == 0:
            self._swanlab.init(
                project=project,
                workspace=workspace,
                experiment_name=experiment_name,
                description=description,
                logdir=self.log_dir,
                mode=mode,
                **kwargs,
            )

    def log_config(self, config: DictConfig) -> None:
        """Saves the config locally and also logs the config to W&B. The config is
        stored in the same directory as the checkpoint. You can
        see an example of the logged config to W&B in the following link:
        https://wandb.ai/capecape/torchtune/runs/6053ofw0/files/torchtune_config_j67sb73v.yaml

        Args:
            config (DictConfig): config to log
        """
        if self._swanlab.get_run():
            resolved = OmegaConf.to_container(config, resolve=True)
            self._swanlab.config.update(resolved)

            output_config_fname = Path(
                os.path.join(
                    config.checkpointer.checkpoint_dir,
                    "torchtune_config.yaml",
                )
            )

            OmegaConf.save(config, output_config_fname)

    def log(self, name: str, data: Scalar, step: int) -> None:
        if self._swanlab.get_run():
            self._swanlab.log({name: data, "global_step": step})

    def log_dict(self, payload: Mapping[str, Scalar], step: int) -> None:
        if self._swanlab.get_run():
            self._swanlab.log({**payload}, step=step)

    def __del__(self) -> None:
        if self._swanlab.get_run():
            self._swanlab.finish()

    def close(self) -> None:
        if self._swanlab.get_run():
            self._swanlab.finish()
