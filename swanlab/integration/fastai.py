from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional

import swanlab
import swanlab.vendor
from swanlab import Callback

_FastaiCallback = swanlab.vendor.fastai.learner.Callback
_fastai_hook = swanlab.vendor.fastai.callback.hook
_fastcore_basics = swanlab.vendor.fastcore.basics


class SwanLabCallback(Callback, _FastaiCallback):
    """
    FastAI callback implementing both ``swanlab.Callback`` and fastai's ``Callback``.

    Usage (recommended):

        import swanlab
        from swanlab.integration.fastai import SwanLabCallback

        cb = SwanLabCallback()
        swanlab.init(project="my-project", callbacks=[cb])
        learn = Learner(..., cbs=cb)
        learn.fit_one_cycle(5)
        swanlab.finish()
    """

    def __init__(
        self,
        *,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        log_dir: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[str] = None,
        tags: Optional[List[str]] = None,
        config: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        if logdir is not None:
            warnings.warn(
                "The `logdir` parameter is deprecated, use `log_dir` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            log_dir = logdir

        Callback.__init__(self)
        _FastaiCallback.__init__(self)

        tags = list(tags) if tags else []
        if "fastai" not in tags:
            tags.append("fastai")

        self._init_kwargs: Dict[str, Any] = {}
        for key, value in [
            ("project", project),
            ("workspace", workspace),
            ("experiment_name", experiment_name),
            ("description", description),
            ("log_dir", log_dir),
            ("mode", mode),
            ("tags", tags),
        ]:
            if value is not None:
                self._init_kwargs[key] = value
        self._init_kwargs.update(kwargs)
        self._init_kwargs.pop("callbacks", None)
        self._init_kwargs.pop("config", None)

        self._pending_config: Dict[str, Any] = dict(config) if config else {}
        self._initialized = False
        self._swanlab_step = 0

    @property
    def name(self) -> str:
        return "swanlab-integration-fastai"

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir: str, path: str, **kwargs: Any) -> None:
        run = self._get_active_run()
        if run is None:
            return
        self._initialized = True
        run.config["FRAMEWORK"] = "fastai"
        self._flush_pending_config(run)

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs: Any) -> None:
        self._initialized = False
        self._pending_config.clear()

    # --- fastai callback hooks ---

    def before_fit(self) -> None:
        self._ensure_init()
        run = self._get_active_run()
        if run is not None:
            configs_log = self._gather_args()
            formatted_config = self._format_config(configs_log)
            run.config.update(formatted_config)

    def after_batch(self) -> None:
        if not self.training:
            return
        self._swanlab_step += 1
        payload: Dict[str, Any] = {
            "train/loss": self.loss.item(),
            "train/step": self._swanlab_step,
        }
        for i, h in enumerate(self.opt.hypers):
            for k, v in h.items():
                payload[f"train/{k}_{i}"] = v
        swanlab.log(payload, step=self._swanlab_step)

    def after_epoch(self) -> None:
        payload: Dict[str, Any] = {}
        for name, value in zip(self.recorder.metric_names, self.recorder.log):
            if value is not None:
                payload[f"epoch/{name}"] = value
        if payload:
            swanlab.log(payload, step=self.epoch)

    # --- helpers ---

    def update_config(self, config: Dict[str, Any]) -> None:
        run = self._get_active_run()
        if run is None:
            self._pending_config.update(config)
            return
        run.config.update(config)

    def _ensure_init(self) -> None:
        if self._initialized:
            return
        run = self._get_active_run()
        if run is not None:
            self._initialized = True
            run.config["FRAMEWORK"] = "fastai"
            self._flush_pending_config(run)
            return
        init_kwargs = dict(self._init_kwargs)
        if self._pending_config:
            init_kwargs["config"] = dict(self._pending_config)
        swanlab.init(callbacks=[self], **init_kwargs)
        self._pending_config.clear()

    def _gather_args(self) -> Dict[str, Any]:
        cb_args = {type(cb).__name__: getattr(cb, "__stored_args__", True) for cb in self.cbs if cb != self}
        args: Dict[str, Any] = {"Learner": self.learn, **cb_args}
        try:
            n_inp = self.dls.train.n_inp
            args["n_inp"] = n_inp
            xb = self.dls.valid.one_batch()[:n_inp]
            args.update(
                {
                    f"input {n + 1} dim {i + 1}": d
                    for n in range(n_inp)
                    for i, d in enumerate(list(_fastcore_basics.detuplify(xb[n]).shape))
                }
            )
        except Exception:
            pass

        with _fastcore_basics.ignore_exceptions():
            args["batch_size"] = self.dls.bs
            args["batch_per_epoch"] = len(self.dls.train)
            args["model_parameters"] = _fastai_hook.total_params(self.model)[0]
            args["device"] = self.dls.device.type
            args["frozen"] = bool(self.opt.frozen_idx)
            args["frozen_idx"] = self.opt.frozen_idx
            args["dataset/tfms"] = f"{self.dls.dataset.tfms}"
            args["dls/after_item"] = f"{self.dls.after_item}"
            args["dls/before_batch"] = f"{self.dls.before_batch}"
            args["dls/after_batch"] = f"{self.dls.after_batch}"
        return args

    @classmethod
    def _format_config(cls, config: Dict[str, Any]) -> Dict[str, Any]:
        result: Dict[str, Any] = {}
        for key, value in config.items():
            if isinstance(value, dict):
                result[key] = cls._format_config(value)
            else:
                result[key] = cls._format_config_value(value)
        return result

    @classmethod
    def _format_config_value(cls, value: Any) -> Any:
        if isinstance(value, list):
            return [cls._format_config_value(item) for item in value]
        elif hasattr(value, "__stored_args__"):
            return {**cls._format_config(value.__stored_args__), "_name": str(value)}
        return value

    def _flush_pending_config(self, run: Any) -> None:
        if not self._pending_config:
            return
        run.config.update(self._pending_config)
        self._pending_config.clear()

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None


__all__ = ["SwanLabCallback"]
