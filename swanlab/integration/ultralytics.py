"""
Docs: https://docs.swanlab.cn/guide_cloud/integration/integration-ultralytics.html

Usage:
------
from ultralytics import YOLO
import swanlab
from swanlab.integration.ultralytics import add_swanlab_callback

model = YOLO("yolov8n.pt")
cb = add_swanlab_callback(model)

swanlab.init(project="my-project", callbacks=[cb])
model.train(data="coco128.yaml", epochs=3)
swanlab.finish()
------

For distributed training, add to ultralytics source:
  your_env/ultralytics/utils/callbacks/base.py → add_integration_callbacks():
    from swanlab.integration.ultralytics import return_swanlab_callback
    callbacks_list.extend([..., return_swanlab_callback()])
"""

from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING, Any, List, Optional

import swanlab
from swanlab import Callback

if TYPE_CHECKING:
    from pathlib import Path

_processed_plots: dict[str, int] = {}


class SwanLabCallback(Callback):
    def __init__(
        self,
        *,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[str] = None,
        tags: Optional[List[str]] = None,
        **kwargs: Any,
    ) -> None:
        self._init_kwargs: dict[str, Any] = {}
        for key, value in [
            ("project", project),
            ("workspace", workspace),
            ("experiment_name", experiment_name),
            ("description", description),
            ("logdir", logdir),
            ("mode", mode),
            ("tags", tags),
        ]:
            if value is not None:
                self._init_kwargs[key] = value
        self._init_kwargs.update(kwargs)
        tags_val = self._init_kwargs.get("tags", [])
        if "ultralytics" not in tags_val:
            tags_val.append("ultralytics")
        self._init_kwargs["tags"] = tags_val

        self._initialized = False
        self._step_counter: Counter = Counter()
        self._model: Any = None

    @property
    def name(self) -> str:
        return "swanlab-integration-ultralytics"

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir: Path, path: str, **kwargs) -> None:
        run = self._get_active_run()
        if run is not None:
            run.config["FRAMEWORK"] = "ultralytics"
        self._initialized = True

    def on_run_finished(self, state: str, error: Optional[str] = None) -> None:
        self._initialized = False
        self._step_counter.clear()
        _processed_plots.clear()

    # --- Ultralytics callback methods ---

    def on_pretrain_routine_start(self, trainer: Any) -> None:
        self._ensure_init()
        run = self._get_active_run()
        if run is None:
            return
        trainer_args = vars(trainer.args)
        if "project" not in self._init_kwargs and getattr(trainer.args, "project", None):
            pass  # project already set during init
        if "experiment_name" not in self._init_kwargs and getattr(trainer.args, "name", None):
            pass
        run.config.update(trainer_args)

    def on_train_epoch_end(self, trainer: Any) -> None:
        run = self._get_active_run()
        if run is None:
            return
        step = trainer.epoch + 1
        swanlab.log(trainer.label_loss_items(trainer.tloss, prefix="train"), step=step)
        swanlab.log(trainer.lr, step=step)
        if trainer.epoch == 1:
            self._log_plots(trainer.plots, step=step, tag="Plots")

    def on_fit_epoch_end(self, trainer: Any) -> None:
        run = self._get_active_run()
        if run is None:
            return
        step = trainer.epoch + 1
        self._step_counter[step] += 1

        if self._step_counter[step] == 2:
            for key, value in trainer.metrics.items():
                swanlab.log({f"FinalValBestModel/{key}": value}, step=step)
        else:
            swanlab.log(trainer.metrics, step=step)

        self._log_plots(trainer.plots, step=step, tag="Train/Plots")
        self._log_plots(trainer.validator.plots, step=step, tag="Train/ValPlots")

        if trainer.epoch == 0:
            try:
                from ultralytics.utils.torch_utils import model_info_for_loggers

                swanlab.log(model_info_for_loggers(trainer), step=step)
            except ImportError:
                pass

    def on_train_end(self, trainer: Any) -> None:
        step = trainer.epoch + 1
        self._log_plots(trainer.plots, step=step, tag="TrainEnd/Plots")
        self._log_plots(trainer.validator.plots, step=step, tag="TrainEnd/ValPlots")

    # --- helpers ---

    def _log_plots(self, plots: dict, step: int, tag: str) -> None:
        run = self._get_active_run()
        if run is None:
            return
        image_list = []
        for name, params in plots.copy().items():
            timestamp = params["timestamp"]
            if _processed_plots.get(name) != timestamp:
                image_list.append(swanlab.Image(str(name), caption=name.stem))
                _processed_plots[name] = timestamp
        if image_list:
            swanlab.log({tag: image_list}, step=step)

    def _ensure_init(self) -> None:
        if self._initialized:
            return
        if self._get_active_run() is not None:
            self._initialized = True
            return
        swanlab.init(callbacks=[self], **self._init_kwargs)

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None


def add_swanlab_callback(
    model: Any,
    *,
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    experiment_name: Optional[str] = None,
    description: Optional[str] = None,
    logdir: Optional[str] = None,
    mode: Optional[str] = None,
    tags: Optional[List[str]] = None,
    **kwargs: Any,
) -> Any:
    cb = SwanLabCallback(
        project=project,
        workspace=workspace,
        experiment_name=experiment_name,
        description=description,
        logdir=logdir,
        mode=mode,
        tags=tags,
        **kwargs,
    )
    callbacks = {
        "on_pretrain_routine_start": cb.on_pretrain_routine_start,
        "on_train_epoch_end": cb.on_train_epoch_end,
        "on_fit_epoch_end": cb.on_fit_epoch_end,
        "on_train_end": cb.on_train_end,
    }
    for event, callback_fn in callbacks.items():
        model.add_callback(event, callback_fn)
    return model


def return_swanlab_callback(
    *,
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    experiment_name: Optional[str] = None,
    description: Optional[str] = None,
    logdir: Optional[str] = None,
    mode: Optional[str] = None,
    tags: Optional[List[str]] = None,
    **kwargs: Any,
) -> dict:
    cb = SwanLabCallback(
        project=project,
        workspace=workspace,
        experiment_name=experiment_name,
        description=description,
        logdir=logdir,
        mode=mode,
        tags=tags,
        **kwargs,
    )
    return {
        "on_pretrain_routine_start": cb.on_pretrain_routine_start,
        "on_train_epoch_end": cb.on_train_epoch_end,
        "on_fit_epoch_end": cb.on_fit_epoch_end,
        "on_train_end": cb.on_train_end,
    }
