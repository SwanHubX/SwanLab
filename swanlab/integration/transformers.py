from __future__ import annotations

import os
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional

import swanlab
import swanlab.vendor
from swanlab import Callback

if TYPE_CHECKING:
    from pathlib import Path


_SINGLE_VALUE_SCALARS = {
    "train_runtime",
    "train_samples_per_second",
    "train_steps_per_second",
    "train_loss",
    "total_flos",
}


def rewrite_logs(logs: Mapping[str, Any]) -> Dict[str, Any]:
    new_logs = {}
    eval_prefix = "eval_"
    test_prefix = "test_"

    for key, value in logs.items():
        if key.startswith(eval_prefix):
            new_logs[f"eval/{key[len(eval_prefix) :]}"] = value
        elif key.startswith(test_prefix):
            new_logs[f"test/{key[len(test_prefix) :]}"] = value
        else:
            new_logs[f"train/{key}"] = value

    return new_logs


class SwanLabCallback(Callback, swanlab.vendor.transformers.trainer_callback.TrainerCallback):
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
        log_config: bool = True,
        **_kwargs: Any,
    ) -> None:
        self._log_config = log_config
        self._init_kwargs: Dict[str, Any] = {}
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

        self._trainer_initialized = False
        self._pending_config: Dict[str, Any] = {}
        self._eval_logged_steps: set[int] = set()

    @property
    def name(self) -> str:
        return "swanlab-integration-transformers"

    def on_run_initialized(self, run_dir: Path, path: str) -> None:
        self._flush_pending_config()

    def on_run_finished(self, state: str, error: Optional[str] = None) -> None:
        self._trainer_initialized = False
        self._pending_config.clear()
        self._eval_logged_steps.clear()

    def setup(self, args: Any, state: Any, model: Any = None) -> None:
        if self._trainer_initialized or not self._is_world_process_zero(state):
            return

        # Auto-init if no active run (backward compat with legacy scripts)
        if self._get_active_run() is None:
            swanlab.init(callbacks=[self], **self._init_kwargs)

        self._trainer_initialized = True
        if self._log_config:
            self._pending_config.update(self._collect_config(args, state, model))
            self._flush_pending_config()

    def update_config(self, config: Dict[str, Any]) -> None:
        self._pending_config.update(config)
        self._flush_pending_config()

    def on_train_begin(self, args: Any, state: Any, control: Any, model: Any = None, **kwargs: Any) -> None:
        self.setup(args, state, model)

    def on_log(
        self,
        args: Any,
        state: Any,
        control: Any,
        model: Any = None,
        logs: Optional[Mapping[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        self.setup(args, state, model)
        if not logs or not self._can_log(state):
            return

        step = self._get_step(state)
        single_value_logs = {
            f"single_value/{key}": value for key, value in logs.items() if key in _SINGLE_VALUE_SCALARS
        }
        metric_logs = {key: value for key, value in logs.items() if key not in _SINGLE_VALUE_SCALARS}

        if single_value_logs:
            self._log(single_value_logs, step)

        rewritten_logs = rewrite_logs(metric_logs)
        if rewritten_logs:
            rewritten_logs["train/global_step"] = step
            self._log(rewritten_logs, step)

    def on_evaluate(
        self, args: Any, state: Any, control: Any, metrics: Optional[Mapping[str, Any]] = None, **kwargs: Any
    ) -> None:
        self.setup(args, state, kwargs.get("model"))
        if metrics is None:
            metrics = kwargs.get("metrics")
        if not metrics or not self._can_log(state):
            return

        step = self._get_step(state)
        if step in self._eval_logged_steps:
            return
        self._eval_logged_steps.add(step)

        rewritten_metrics = rewrite_logs(metrics)
        if rewritten_metrics:
            self._log(rewritten_metrics, step)

    def on_save(self, args: Any, state: Any, control: Any, **kwargs: Any) -> None:
        if not self._can_log(state):
            return
        run = self._get_active_run()
        if run is None:
            return

        step = self._get_step(state)
        output_dir = getattr(args, "output_dir", None)
        if output_dir:
            run.config["transformers_last_checkpoint"] = os.path.join(output_dir, f"checkpoint-{step}")

    def on_predict(
        self,
        args: Any,
        state: Any,
        control: Any,
        metrics: Optional[Mapping[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        self.setup(args, state, kwargs.get("model"))
        if not metrics or not self._can_log(state):
            return

        self._log(rewrite_logs(metrics), self._get_step(state))

    @staticmethod
    def _is_world_process_zero(state: Any) -> bool:
        return bool(getattr(state, "is_world_process_zero", True))

    @staticmethod
    def _get_step(state: Any) -> int:
        step = getattr(state, "global_step", 0)
        return step if isinstance(step, int) and step >= 0 else 0

    @staticmethod
    def _get_active_run():
        try:
            return swanlab.get_run()
        except RuntimeError:
            return None

    def _can_log(self, state: Any) -> bool:
        return self._is_world_process_zero(state) and self._get_active_run() is not None

    def _collect_config(self, args: Any, state: Any, model: Any = None) -> Dict[str, Any]:
        config: Dict[str, Any] = {"FRAMEWORK": "transformers"}

        model_config = getattr(model, "config", None)
        if isinstance(model_config, dict):
            config.update(model_config)
        else:
            model_config_to_dict = getattr(model_config, "to_dict", None)
            if callable(model_config_to_dict):
                model_config_dict = model_config_to_dict()
                if isinstance(model_config_dict, Mapping):
                    config.update(model_config_dict)

        args_to_dict = getattr(args, "to_dict", None)
        if callable(args_to_dict):
            args_config = args_to_dict()
            if isinstance(args_config, Mapping):
                config.update(args_config)

        peft_config = getattr(model, "peft_config", None)
        if peft_config is not None:
            config["peft_config"] = peft_config

        trial_name = getattr(state, "trial_name", None)
        if trial_name is not None:
            config["transformers_trial_name"] = trial_name

        num_parameters = getattr(model, "num_parameters", None)
        if callable(num_parameters):
            config["model_num_parameters"] = num_parameters()

        trainable_parameters = getattr(model, "get_nb_trainable_parameters", None)
        if callable(trainable_parameters):
            parameter_counts = trainable_parameters()
            if isinstance(parameter_counts, tuple) and len(parameter_counts) == 2:
                trainable_params, all_params = parameter_counts
                config["peft_model_trainable_params"] = trainable_params
                config["peft_model_all_param"] = all_params

        return config

    def _flush_pending_config(self) -> None:
        if not self._pending_config:
            return
        run = self._get_active_run()
        if run is None:
            return
        run.config.update(self._pending_config)
        self._pending_config.clear()

    def _log(self, payload: Mapping[str, Any], step: int) -> None:
        swanlab.log(dict(payload), step=step)
