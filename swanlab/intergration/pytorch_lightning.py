import os
import importlib.util
from typing import Any, Dict, Optional, Union, Mapping
from argparse import Namespace

import packaging.version

if importlib.util.find_spec("lightning"):
    import lightning.pytorch as pl

    from lightning.pytorch.loggers.logger import Logger, rank_zero_experiment

    from lightning.pytorch.utilities import rank_zero_only
elif importlib.util.find_spec("pytorch_lightning"):
    import pytorch_lightning as pl

    if packaging.version.parse(pl.__version__) < packaging.version.parse("1.7"):
        from pytorch_lightning.loggers.base import (
            LightningLoggerBase as Logger,
            rank_zero_experiment,
        )
    else:
        from pytorch_lightning.loggers.logger import (
            Logger,
            rank_zero_experiment,
        )

    from pytorch_lightning.utilities import rank_zero_only
else:
    raise RuntimeError(
        "This contrib module requires PyTorch Lightning to be installed. "
        "Please install it with command: \n pip install pytorch-lightning"
        "or \n pip install lightning"
    )

from ..data.run import SwanLabRun
from lightning_fabric.utilities.logger import _add_prefix, _convert_params, _sanitize_callable_params


"""
swanlab_logger = SwanLabLogger(
    project="example_project",
    description="example_description",
)
"""


class SwanLabLogger(Logger):

    LOGGER_JOIN_CHAR = "-"

    def __init__(
        self,
        project: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        prefix: str = "",
        logdir: Optional[str] = None,
        cloud: Optional[bool] = None,
        **kwargs: Any,
    ):
        super().__init__()
        self._prefix = prefix
        self._experiment = None

        self._swanlab_init: Dict[str, Any] = {
            "project": project,
            "experiment_name": experiment_name,
            "description": description,
            "logdir": logdir,
            "cloud": cloud,
        }

        self._swanlab_init.update(**kwargs)

        self._project = self._swanlab_init.get("project")
        self._experiment_name = self._swanlab_init.get("experiment_name")
        self._description = self._swanlab_init.get("decsription")
        self._save_dir = self._swanlab_init.get("logdir")
        self._cloud = self._swanlab_init.get("cloud")

    @property
    @rank_zero_experiment
    def experiment(self) -> SwanLabRun:
        """创建实验"""
        import swanlab

        if self._experiment is None:
            self._experiment = swanlab.init(**self._swanlab_init)

        return self._experiment

    @rank_zero_only
    def log_hyperparams(self, params: Union[Dict[str, Any], Namespace]) -> None:
        """记录超参数(config)"""
        params = _convert_params(params)
        params = _sanitize_callable_params(params)

        config: dict = self.experiment.config

        config.update(params, allow_val_change=True)

    @rank_zero_only
    def log_metrics(self, metrics: Mapping[str, float], step: Optional[int] = None) -> None:
        """记录标量指标"""

        assert rank_zero_only.rank == 0, "experiment tried to log from global_rank != 0"

        metrics = _add_prefix(metrics, self._prefix, self.LOGGER_JOIN_CHAR)
        if step is not None:
            self.experiment.log(metrics, step=int(step))
        else:
            self.experiment.log(metrics)

    @property
    def save_dir(self) -> Optional[str]:
        return self._save_dir

    @property
    def name(self) -> Optional[str]:
        return self._project

    @rank_zero_only
    def finalize(self, status: str) -> None:
        """结束记录"""
        # import swanlab

        # if status != "success":
        #     # Currently, checkpoints only get logged on success
        #     return

        # if self._experiment is not None:
        #     swanlab.finish()
        pass

    @property
    def version(self) -> Optional[str]:
        # don't create an experiment if we don't have one
        return self.experiment.settings.run_id if self._experiment is not None else None
