from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Union

import swanlab
import swanlab.vendor
from swanlab import Callback

LogStrategy = Literal["epoch", "batch"]

_KerasCallback = swanlab.vendor.keras.callbacks.Callback


class SwanLabCallback(Callback, _KerasCallback):
    def __init__(
        self,
        log_freq: Union[LogStrategy, int] = "epoch",
        initial_global_step: int = 0,
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
            import warnings

            warnings.warn(
                "The `logdir` parameter is deprecated, use `log_dir` instead.", DeprecationWarning, stacklevel=2
            )
            log_dir = logdir
        _KerasCallback.__init__(self)

        if log_freq == "batch":
            log_freq = 1
        if isinstance(log_freq, int) and log_freq <= 0:
            raise ValueError("log_freq must be a positive integer, 'batch', or 'epoch'")

        self.logging_batch_wise = isinstance(log_freq, int)
        self.log_freq: Optional[int] = log_freq if isinstance(log_freq, int) else None
        self.global_batch = 0
        self.global_step = initial_global_step
        self._initialized = False
        self._pending_config: Dict[str, Any] = dict(config) if config else {}
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

    @property
    def name(self) -> str:
        return "swanlab-integration-keras"

    # --- swanlab.Callback hooks ---

    def on_run_initialized(self, run_dir, path, **kwargs) -> None:
        run = self._get_active_run()
        if run is None:
            return
        self._initialized = True
        run.config["FRAMEWORK"] = "keras"
        # _pending_config already passed via swanlab.init(config=...) in _ensure_init,
        # so no need to flush again here.

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs) -> None:
        self._initialized = False
        self._pending_config.clear()

    # --- Keras callback hooks ---

    def update_config(self, config: Dict[str, Any]) -> None:
        run = self._get_active_run()
        if run is None:
            self._pending_config.update(config)
            return
        run.config.update(config)

    def on_epoch_end(self, epoch: int, logs: Optional[Dict[str, Any]] = None) -> None:
        payload = {f"epoch/{key}": value for key, value in (logs or {}).items()}
        payload["epoch/epoch"] = epoch

        learning_rate = self._get_lr()
        if learning_rate is not None:
            payload["epoch/learning_rate"] = learning_rate

        self._log(payload, step=epoch)

    def on_batch_end(self, batch: int, logs: Optional[Dict[str, Any]] = None) -> None:
        self.global_step += 1
        if not self.logging_batch_wise or self.log_freq is None or batch % self.log_freq != 0:
            return

        payload = {f"batch/{key}": value for key, value in (logs or {}).items()}
        payload["batch/batch_step"] = self.global_batch

        learning_rate = self._get_lr()
        if learning_rate is not None:
            payload["batch/learning_rate"] = learning_rate

        self._log(payload, step=self.global_step)
        self.global_batch += self.log_freq

    def on_train_batch_end(self, batch: int, logs: Optional[Dict[str, Any]] = None) -> None:
        self.on_batch_end(batch, logs)

    # --- helpers ---

    def _ensure_init(self) -> None:
        if self._initialized:
            return
        run = self._get_active_run()
        if run is not None:
            self._initialized = True
            run.config["FRAMEWORK"] = "keras"
            self._flush_pending_config(run)
            return
        init_kwargs = dict(self._init_kwargs)
        if self._pending_config:
            init_kwargs["config"] = dict(self._pending_config)
        swanlab.init(callbacks=[self], **init_kwargs)
        self._pending_config.clear()

    def _log(self, payload: Dict[str, Any], step: int) -> None:
        self._ensure_init()
        if self._get_active_run() is None:
            return
        swanlab.log(payload, step=step)

    def _flush_pending_config(self, run: Any) -> None:
        if not self._pending_config:
            return
        run.config.update(self._pending_config)
        self._pending_config.clear()

    def _get_lr(self) -> Optional[float]:
        model = getattr(self, "model", None)
        optimizer = getattr(model, "optimizer", None)
        if optimizer is None:
            return None

        learning_rate = getattr(optimizer, "learning_rate", None)
        if learning_rate is None:
            learning_rate = getattr(optimizer, "lr", None)
        if learning_rate is None:
            return None

        try:
            value = learning_rate(self.global_step) if callable(learning_rate) else learning_rate
            return self._to_float(value)
        except Exception:
            return None

    @staticmethod
    def _to_float(value: Any) -> float:
        numpy = getattr(value, "numpy", None)
        if callable(numpy):
            value = numpy()
        item = getattr(value, "item", None)
        if callable(item):
            value = item()
        return float(value)

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None


SwanLabLogger = SwanLabCallback

__all__ = ["LogStrategy", "SwanLabCallback", "SwanLabLogger"]
