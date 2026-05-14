from __future__ import annotations

from typing import Any, Dict, List, Optional

import swanlab
import swanlab.vendor
from swanlab import Callback

_Sb3BaseCallback = swanlab.vendor.stable_baselines3.common.callbacks.BaseCallback  # type: ignore
_Sb3KVWriter = swanlab.vendor.stable_baselines3.common.logger.KVWriter  # type: ignore
_Sb3Logger = swanlab.vendor.stable_baselines3.common.logger.Logger  # type: ignore


class _SwanLabOutputFormat(_Sb3KVWriter):
    """Intercept SB3 logger writes and forward scalar metrics to swanlab."""

    def write(
        self,
        key_values: Dict[str, Any],
        key_excluded: Dict[str, Any],
        step: int = 0,
    ) -> None:
        logs: Dict[str, Any] = {}
        for key, value in key_values.items():
            if isinstance(value, (int, float)):
                logs[key] = value
        if logs:
            swanlab.log(logs, step=step)


class SwanLabCallback(Callback, _Sb3BaseCallback):
    """
    Stable Baselines3 callback implementing both ``swanlab.Callback``
    and SB3's ``BaseCallback``.

    Usage (recommended):

        import swanlab
        from swanlab.integration.sb3 import SwanLabCallback

        cb = SwanLabCallback()
        swanlab.init(project="my-project", callbacks=[cb])
        model = PPO("MlpPolicy", env, verbose=1)
        model.learn(total_timesteps=25000, callback=cb)
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
        verbose: int = 0,
        **kwargs: Any,
    ) -> None:
        if logdir is not None:
            import warnings

            warnings.warn(
                "The `logdir` parameter is deprecated, use `log_dir` instead.", DeprecationWarning, stacklevel=2
            )
            log_dir = logdir

        Callback.__init__(self)
        _Sb3BaseCallback.__init__(self, verbose=verbose)

        tags = list(tags) if tags else []
        if "stable_baselines3" not in tags:
            tags.append("stable_baselines3")

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

        self._initialized = False
        self._pending_config: Dict[str, Any] = {}

    @property
    def name(self) -> str:
        return "swanlab-integration-sb3"

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir: str, path: str, **kwargs: Any) -> None:
        self._initialized = True
        run = self._get_active_run()
        if run is not None:
            run.config["FRAMEWORK"] = "stable_baselines3"
            self._flush_pending_config()

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs: Any) -> None:
        self._initialized = False
        self._pending_config.clear()

    # --- SB3 callback hooks ---

    def _on_training_start(self) -> None:
        self._ensure_init()
        self._collect_config()

        # Install SwanLab output format into SB3 logger.
        # Preserve existing output formats so console printing is not lost.
        existing_logger = getattr(self.model, "logger", None)
        output_formats: List[Any] = [_SwanLabOutputFormat()]
        if existing_logger is not None and hasattr(existing_logger, "output_formats"):
            output_formats.extend(existing_logger.output_formats)

        new_logger = _Sb3Logger(folder=None, output_formats=output_formats)
        self.model.set_logger(new_logger)

    def _on_step(self) -> bool:
        return True

    # --- helpers ---

    def _ensure_init(self) -> None:
        if self._initialized:
            return
        run = self._get_active_run()
        if run is not None:
            self._initialized = True
            run.config["FRAMEWORK"] = "stable_baselines3"
            self._flush_pending_config()
            return
        init_kwargs = dict(self._init_kwargs)
        if self._pending_config:
            init_kwargs["config"] = dict(self._pending_config)
        swanlab.init(callbacks=[self], **init_kwargs)
        self._pending_config.clear()

    def update_config(self, config: Dict[str, Any]) -> None:
        """Update swanlab config (backward-compatible helper)."""
        run = self._get_active_run()
        if run is None:
            self._pending_config.update(config)
            return
        run.config.update(config)

    def _collect_config(self) -> None:
        """Collect model and training hyperparameters into config."""
        run = self._get_active_run()
        if run is None:
            return

        config: Dict[str, Any] = {}
        if hasattr(self, "model") and self.model is not None:
            model = self.model
            config["algo"] = type(model).__name__
            for key, value in model.__dict__.items():
                if key.startswith("_"):
                    continue
                if isinstance(value, (float, int, str, bool)):
                    config[key] = value
                elif hasattr(value, "__name__"):
                    config[key] = value.__name__

            policy = getattr(model, "policy", None)
            if policy is not None and hasattr(policy, "net_arch"):
                config["net_arch"] = str(policy.net_arch)

        run.config.update(config)

    def _flush_pending_config(self) -> None:
        if not self._pending_config:
            return
        run = self._get_active_run()
        if run is None:
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
