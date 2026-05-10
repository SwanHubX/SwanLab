"""
Docs: https://docs.swanlab.cn/guide_cloud/integration/integration-mmengine.html

Usage:
------
# In your mmengine config:
custom_imports = dict(
    imports=["swanlab.integration.mmengine"], allow_failed_imports=False
)
vis_backends = [
    dict(
        type="SwanLabVisBackend",
        init_kwargs={
            "project": "YourProject",
            "experiment_name": "YourExperiment",
        },
    ),
]
------
"""

from __future__ import annotations

import os
from typing import Any, Optional, Sequence, Union

import swanlab
import swanlab.vendor
from swanlab import Callback

_BaseVisBackend = swanlab.vendor.mmengine.visualization.vis_backend.BaseVisBackend
_VISBACKENDS = swanlab.vendor.mmengine.registry.VISBACKENDS


def _repack_dict(d: dict, prefix: str = "") -> dict:
    result: dict = {}
    for key, value in d.items():
        key = str(key)
        full_key = f"{prefix}/{key}" if prefix else key
        if isinstance(value, dict):
            result.update(_repack_dict(value, full_key))
        elif isinstance(value, (list, tuple)):
            if all(not isinstance(item, dict) for item in value):
                result[full_key] = value
            else:
                for i, item in enumerate(value):
                    result.update(_repack_dict(item, f"{full_key}[{i}]"))
        else:
            result[full_key] = value
    return result


@_VISBACKENDS.register_module()
class SwanlabVisBackend(Callback, _BaseVisBackend):
    """SwanLab visualization backend for mmengine.

    Args:
        save_dir: Root directory for SwanLab files.
        init_kwargs: Keyword arguments passed to ``swanlab.init()``.
    """

    def __init__(
        self,
        save_dir: Optional[str] = None,
        init_kwargs: Optional[dict] = None,
    ) -> None:
        self._save_dir = save_dir
        self._env_initialized = False
        self._init_kwargs = init_kwargs or {}
        self._swanlab_initialized = False

    @property
    def name(self) -> str:
        return "swanlab-integration-mmengine"

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir, path, **kwargs) -> None:
        run = self._get_active_run()
        if run is not None:
            run.config["FRAMEWORK"] = "mmengine"
        self._swanlab_initialized = True

    def on_run_finished(self, state: str, error: Optional[str] = None) -> None:
        self._swanlab_initialized = False

    # --- mmengine BaseVisBackend interface ---

    def _init_env(self) -> None:
        if self._save_dir is not None:
            os.makedirs(self._save_dir, exist_ok=True)
        self._init_kwargs.setdefault("log_dir", self._save_dir)

        if self._get_active_run() is not None:
            self._swanlab_initialized = True
        else:
            swanlab.init(callbacks=[self], **self._init_kwargs)

        self._env_initialized = True

    @property
    def experiment(self) -> Any:
        if not self._env_initialized:
            self._init_env()
        return swanlab

    def add_config(self, config: Any, **kwargs) -> None:
        if not self._env_initialized:
            self._init_env()
        config_dict = config.to_dict() if hasattr(config, "to_dict") else dict(config)
        run = self._get_active_run()
        if run is not None:
            run.config.update(_repack_dict(config_dict))

    def add_graph(self, model: Any, data_batch: Sequence[dict], **kwargs) -> None:
        pass

    def add_image(self, name: str, image: Any, step: int = 0, **kwargs) -> None:
        if not self._env_initialized:
            self._init_env()
        swanlab.log({name: swanlab.Image(image)}, step=step)

    def add_scalar(self, name: str, value: Union[int, float], step: int = 0, **kwargs) -> None:
        if not self._env_initialized:
            self._init_env()
        swanlab.log({name: value}, step=step)

    def add_scalars(
        self,
        scalar_dict: dict,
        step: int = 0,
        file_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        if not self._env_initialized:
            self._init_env()
        swanlab.log(scalar_dict, step=step)

    def close(self) -> None:
        pass

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None
