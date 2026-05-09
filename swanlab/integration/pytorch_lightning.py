from __future__ import annotations

import os
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Union

import swanlab
import swanlab.vendor
from swanlab import Callback
from swanlab.sdk.internal.pkg import console

_LightningLogger = swanlab.vendor.lightning.pytorch.loggers.Logger
_rank_zero_only = swanlab.vendor.lightning.pytorch.utilities.rank_zero_only
_rank_zero_warn = getattr(swanlab.vendor.lightning.pytorch.utilities, "rank_zero_warn", console.warning)


class SwanLabLogger(Callback, _LightningLogger):
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[str] = None,
        save_dir: Optional[Union[str, Path]] = ".",
        tags: Optional[List[str]] = None,
        id: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        _LightningLogger.__init__(self)

        self._experiment = None
        self._initialized = False
        self._pending_config: Dict[str, Any] = {}
        self._project = project
        self._experiment_name = experiment_name
        self._save_dir = os.fspath(save_dir) if save_dir is not None else None
        self._init_kwargs: Dict[str, Any] = {}

        for key, value in [
            ("project", project),
            ("workspace", workspace),
            ("experiment_name", experiment_name),
            ("description", description),
            ("logdir", logdir),
            ("mode", mode),
            ("tags", tags),
            ("id", id),
        ]:
            if value is not None:
                self._init_kwargs[key] = value
        self._init_kwargs.update(kwargs)
        self._init_kwargs.pop("callbacks", None)

    @property
    def name(self) -> str:
        return "swanlab-integration-pytorch-lightning"

    @property
    def version(self) -> Optional[str]:
        run = self._get_active_run()
        if run is not None:
            return getattr(run, "id", None)
        if self._experiment is not None:
            return getattr(self._experiment, "id", None)
        return None

    @property
    def save_dir(self) -> Optional[str]:
        return self._save_dir

    @property
    def log_dir(self) -> Optional[str]:
        return self._save_dir

    @property
    def experiment(self):
        self._ensure_init()
        run = self._get_active_run()
        if run is not None:
            self._experiment = run
            return run
        return self

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir, path, **kwargs) -> None:
        run = self._get_active_run()
        if run is None:
            return
        self._initialized = True
        self._experiment = run
        run.config["FRAMEWORK"] = "pytorch_lightning"
        self._flush_pending_config(run)

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs) -> None:
        self._initialized = False
        self._experiment = None
        self._pending_config.clear()

    # --- Lightning Logger interface ---

    def update_config(self, config: Dict[str, Any]) -> None:
        run = self._get_active_run()
        if run is None:
            self._pending_config.update(config)
            return
        run.config.update(config)

    @_rank_zero_only
    def log_hyperparams(self, params: Union[Mapping[str, Any], Namespace], *args: Any, **kwargs: Any) -> None:
        self._ensure_init()
        run = self._get_active_run()
        if run is None:
            return
        run.config.update(self._params_to_dict(params))

    @_rank_zero_only
    def log_metrics(self, metrics: Mapping[str, Any], step: Optional[int] = None) -> None:
        self._ensure_init()
        if self._get_active_run() is None:
            return
        swanlab.log(dict(metrics), step=step)

    @_rank_zero_only
    def log_image(self, key: str, images: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        self._log_media_list(swanlab.Image, key, images, step, **kwargs)

    @_rank_zero_only
    def log_audio(self, key: str, audios: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        self._log_media_list(swanlab.Audio, key, audios, step, **kwargs)

    @_rank_zero_only
    def log_text(self, key: str, texts: List[Any], step: Optional[int] = None, **kwargs: Any) -> None:
        self._log_media_list(swanlab.Text, key, texts, step, **kwargs)

    @_rank_zero_only
    def save(self) -> None:
        pass

    @_rank_zero_only
    def finalize(self, status: Optional[str] = None) -> None:
        if status is not None and status != "success":
            _rank_zero_warn(f"SwanLabLogger received non-success Lightning finalize status: {status}")

    # --- helpers ---

    def _ensure_init(self) -> None:
        if self._initialized:
            return
        run = self._get_active_run()
        if run is not None:
            self._initialized = True
            self._experiment = run
            run.config["FRAMEWORK"] = "pytorch_lightning"
            self._flush_pending_config(run)
            return
        init_kwargs = dict(self._init_kwargs)
        if self._pending_config:
            init_kwargs["config"] = dict(self._pending_config)
        swanlab.init(callbacks=[self], **init_kwargs)
        self._pending_config.clear()

    def _flush_pending_config(self, run: Any) -> None:
        if not self._pending_config:
            return
        run.config.update(self._pending_config)
        self._pending_config.clear()

    def _log_media_list(self, factory: Any, key: str, values: List[Any], step: Optional[int], **kwargs: Any) -> None:
        if not isinstance(values, list):
            raise TypeError(f'Expected a list for "{key}", found {type(values).__name__}')

        value_count = len(values)
        for arg_key, arg_values in kwargs.items():
            if len(arg_values) != value_count:
                raise ValueError(f"Expected {value_count} items but only found {len(arg_values)} for {arg_key}")

        per_value_kwargs = [
            {
                arg_key: arg_values[index]
                if isinstance(arg_values, (list, tuple)) and len(arg_values) == value_count
                else arg_values
                for arg_key, arg_values in kwargs.items()
            }
            for index in range(value_count)
        ]
        self.log_metrics(
            {key: [factory(value, **value_kwargs) for value, value_kwargs in zip(values, per_value_kwargs)]}, step
        )

    @staticmethod
    def _params_to_dict(params: Union[Mapping[str, Any], Namespace]) -> Dict[str, Any]:
        if isinstance(params, Namespace):
            params = vars(params)
        elif hasattr(params, "__dict__"):
            params = vars(params)
        elif not isinstance(params, Mapping):
            params = {}

        sanitized = {}
        for key, value in dict(params).items():
            if callable(value):
                value = getattr(value, "__name__", str(value))
            sanitized[key] = value
        return sanitized

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None


__all__ = ["SwanLabLogger"]
