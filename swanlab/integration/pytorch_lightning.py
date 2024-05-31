"""
Docs:https://docs.swanlab.cn/zh/guide_cloud/integration/integration-pytorch-lightning.html
"""

import os
import importlib.util
from typing import Any, Dict, Optional, Union, Mapping, List
from argparse import Namespace
import packaging.version
from pathlib import Path

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
import swanlab
from ..data.run import SwanLabRun
from lightning_fabric.utilities.logger import _add_prefix, _convert_params, _sanitize_callable_params


"""
swanlab_logger = SwanLabLogger(
    project="example_project",
    description="example_description",
)
"""


class SwanLabLogger(Logger):
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[str] = None,
        save_dir: Union[str, Path] = ".",
        **kwargs: Any,
    ):
        super().__init__()

        self._swanlab_init: Dict[str, Any] = {
            "project": project,
            "workspace": workspace,
            "experiment_name": experiment_name,
            "description": description,
            "logdir": logdir,
            "mode": mode,
        }

        self._swanlab_init.update(**kwargs)

        self._project = self._swanlab_init.get("project")
        self._workspace = self._swanlab_init.get("workspace")
        self._experiment_name = self._swanlab_init.get("experiment_name")
        self._description = self._swanlab_init.get("decsription")
        self._logdir = self._swanlab_init.get("logdir")
        self._mode = self._swanlab_init.get("mode")

        if save_dir is not None:
            save_dir = os.fspath(save_dir)
        self._save_dir = save_dir

    @property
    @rank_zero_experiment
    def experiment(self) -> SwanLabRun:
        """创建实验"""
        if swanlab.get_run() is None:
            self._experiment = swanlab.init(**self._swanlab_init)
        else:
            self._experiment = swanlab.get_run()

        return self._experiment

    @rank_zero_only
    def log_hyperparams(self, params: Union[Dict[str, Any], Namespace]) -> None:
        """记录超参数(config)"""
        params = _convert_params(params)
        params = _sanitize_callable_params(params)

        self.experiment.config.update(params)

    @rank_zero_only
    def log_metrics(self, metrics: Mapping[str, float], step: Optional[int] = None) -> None:
        """记录标量指标"""

        assert rank_zero_only.rank == 0, "experiment tried to log from global_rank != 0"

        if step is not None:
            self.experiment.log(metrics, step=step)
        else:
            self.experiment.log(metrics)

    @rank_zero_only
    def log_image(self, key: str, images: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        """Log images (tensors, numpy arrays, PIL Images or file paths).

        Optional kwargs are lists passed to each image (ex: caption).

        """
        if not isinstance(images, list):
            raise TypeError(f'Expected a list as "images", found {type(images)}')
        n = len(images)
        for k, v in kwargs.items():
            if len(v) != n:
                raise ValueError(f"Expected {n} items but only found {len(v)} for {k}")
        kwarg_list = [{k: kwargs[k][i] for k in kwargs} for i in range(n)]

        import swanlab

        metrics = {key: [swanlab.Image(img, **kwarg) for img, kwarg in zip(images, kwarg_list)]}
        self.log_metrics(metrics, step)  # type: ignore[arg-type]

    @rank_zero_only
    def log_audio(self, key: str, audios: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        r"""Log audios (numpy arrays, or file paths).

        Args:
            key: The key to be used for logging the audio files
            audios: The list of audio file paths, or numpy arrays to be logged
            step: The step number to be used for logging the audio files
            \**kwargs: Optional kwargs are lists passed to each ``swanlab.Audio`` instance (ex: caption, sample_rate).

        Optional kwargs are lists passed to each audio (ex: caption, sample_rate).

        """
        if not isinstance(audios, list):
            raise TypeError(f'Expected a list as "audios", found {type(audios)}')
        n = len(audios)
        for k, v in kwargs.items():
            if len(v) != n:
                raise ValueError(f"Expected {n} items but only found {len(v)} for {k}")
        kwarg_list = [{k: kwargs[k][i] for k in kwargs} for i in range(n)]

        import swanlab

        metrics = {key: [swanlab.Audio(audio, **kwarg) for audio, kwarg in zip(audios, kwarg_list)]}
        self.log_metrics(metrics, step)  # type: ignore[arg-type]

    @rank_zero_only
    def log_text(self, key: str, texts: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        r"""Log texts (numpy arrays, or file paths).

        Args:
            key: The key to be used for logging the string
            audios: The list of string to be logged
            step: The step number to be used for logging the string
            \**kwargs: Optional kwargs are lists passed to each ``swanlab.Audio`` instance (ex: caption, sample_rate).

        Optional kwargs are lists passed to each text (ex: caption).

        """
        if not isinstance(texts, list):
            raise TypeError(f'Expected a list as "texts", found {type(texts)}')
        n = len(texts)
        for k, v in kwargs.items():
            if len(v) != n:
                raise ValueError(f"Expected {n} items but only found {len(v)} for {k}")
        kwarg_list = [{k: kwargs[k][i] for k in kwargs} for i in range(n)]

        import swanlab

        metrics = {key: [swanlab.Text(text, **kwarg) for text, kwarg in zip(texts, kwarg_list)]}
        self.log_metrics(metrics, step)  # type: ignore[arg-type]

    @property
    def save_dir(self) -> Optional[str]:
        return self._save_dir

    @property
    def name(self) -> Optional[str]:
        return self._project

    @rank_zero_only
    def finalize(self, status: str) -> None:
        if status != "success":
            return

    @property
    def version(self) -> Optional[str]:
        # don't create an experiment if we don't have one
        return self.experiment.settings.run_id if self._experiment is not None else None
